#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#devtools::install_github('betsig/ipath')
library(shiny)
library(rmarkdown)
library(viridis)
library(viridisLite)

options(stringsAsFactors = FALSE)
options(shiny.maxRequestSize=100*1024^2)

#### ipath source code .... ####

library(stringr);library(RColorBrewer);library(httr);library(rsvg);library(reshape2);library(plyr);library(ggplot2);library(ggforce)


#' Write an iPath3 metabolic pathways map to pdf
#' @param ipath_data data.frame formatted for sending to ipath
#' @param file_name pdf file name to export plot.
#' @param default_color Hex color value. Non-selected parts of the map will use this color by default.
#' @param default_width Default path/edge width (px) for the non-selected parts of the map.
#' @param default_radius Compound radius (px) for the non-selected parts of the map.
#' @param textbox_textcol hex value for text (in box) color
#' @param textbox_col hex value for box color
#' @export
#' @import httr
#' @importFrom rsvg rsvg_pdf
#' @author Beth Signal
#' @examples
print_ipath_pdf = function(ipath_data, file_name=NULL, default_width = 2,
                           default_color = '#bdbdbd',default_radius = 5,
                           textbox_textcol = NULL, textbox_col = NULL){
  
  params = to_parameters(ipath_data,
                         default_width = default_width,
                         default_color = default_color,
                         default_radius = default_radius)
  
  ipath_result = httr::POST("https://pathways.embl.de/mapping.cgi?map=metabolic",
                            body = params, encode = "form")
  
  pathway_map = httr::content(ipath_result, as="text", encoding = "UTF-8")
  pathway_map = colour_on_top(pathway_map, ipath_data)
  
  if(!is.null(textbox_textcol) | !is.null(textbox_col)){
    pathway_map = recolor_textbox(pathway_map, box_color = textbox_col, text_color = textbox_textcol)
  }
  
  #rsvg::rsvg_pdf(charToRaw(pathway_map), file = file_name)
  # scale the svg so it's easily viewble
  pathway_map = gsub("height='2250' width='3774' viewBox=\"-10 -10 3794 2270\"", "height='16cm' width='26cm' viewBox=\"0 0 26.0cm 16.0cm\"", pathway_map)
  pathway_map = gsub("<style type='text/css'>", "<g transform=\"translate(0, 0) scale(0.265) \"><style type='text/css'>", pathway_map)
  pathway_map = gsub("</g>\n</svg>", "</g>\n</g>\n</svg>", pathway_map)
  
  # remove crap
  #pathway_map = gsub("[\n][<]text(.*?)Kanehisa(.*?)text[>]", "", pathway_map)
  #pathway_map = gsub("[\n][<]text(.*?)5/11/17(.*?)text[>]", "", pathway_map)
  
  all_newlines = str_locate_all(pathway_map, "\n")
  
  rm_loc = str_locate_all(pathway_map, "5/11/17")
  diffs = all_newlines[[1]][,1] - rm_loc[[1]][1]
  str_sub(pathway_map, all_newlines[[1]][which(diffs == max(diffs[diffs < 0])),1],all_newlines[[1]][which(diffs == min(diffs[diffs > 0])),1]-1) = ""
  
  rm_loc = str_locate_all(pathway_map, "Kanehisa")
  diffs = all_newlines[[1]][,1] - rm_loc[[1]][1]
  str_sub(pathway_map, all_newlines[[1]][which(diffs == max(diffs[diffs < 0])),1],all_newlines[[1]][which(diffs == min(diffs[diffs > 0])),1]-1) = ""
  
  if(!is.null(file_name)){
    write.table(pathway_map, file = file_name, row.names = F, col.names = F, quote = F, sep="\n")
    #svg_xml = data.table::fread(file_name,  data.table=FALSE, sep='\n', header = F)[,1]
  }
  
  svg_xml = unlist(str_split(pathway_map, "\n"))
  
  return(svg_xml)
  
}

#' Convert svg data to tables for ggplot
#' @param svg_xml data.frame containing svg data
#' @export
#' @import stringr
#' @author Beth Signal
#' @examples
svg_to_table = function(svg_xml){
  
  #svg_xml = data.table::fread(svg_file, data.table = FALSE, header = FALSE)
  
  layer = c(1:5)
  layer_loc = vector()
  for(i in layer){
    layer_loc[i] = grep(paste0("layer", i), svg_xml)
  }
  
  layer_locations = data.frame(layer,layer_loc)
  
  ## Layer 1: paths ##
  layer_1_paths = grep("path", svg_xml[layer_locations$layer_loc[1]:layer_locations$layer_loc[2]], value = TRUE)
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_1_paths, "[ ]+"), "[[",3))))
  stroke = gsub("'", "",gsub('"',"",(gsub('stroke=',"",unlist(lapply(stringr::str_split(layer_1_paths, "[ ]+"), "[[",5))))))
  stroke_width = as.numeric(gsub("'", "",gsub('"',"",(gsub('stroke-width=',"",unlist(lapply(stringr::str_split(layer_1_paths, "[ ]+"), "[[",7)))))))
  
  path_locs = gsub("'", "",gsub('"',"",gsub('/>',"",gsub('d=',"",unlist(lapply(stringr::str_split(lapply(stringr::str_split(layer_1_paths, "d="), "[[",2), ">"), "[[",1))))))
  
  layer_1 = data.frame(opacity, stroke, stroke_width, id = path_locs)
  
  path_locs_CZ = path_locs[grepl("C", path_locs) | grepl("Z", path_locs)]
  path_locs_L = path_locs[!(grepl("C", path_locs) | grepl("Z", path_locs))]
  
  coords_CZ = do.call("rbind", lapply(unique(path_locs_CZ), function(x) cbind(id = x, make_svg_df(x))))
  layer_1_curves = layer_1[layer_1$id %in% path_locs_CZ,]
  m = match(coords_CZ$id, layer_1_curves$id)
  layer_1_curves = cbind(layer_1_curves[m,], coords_CZ[,-which(colnames(coords_CZ) == "id")])
  color_order = as.data.frame(table(layer_1_curves$stroke))
  base_col = color_order$Var1[which.max(color_order$Freq)]
  layer_1_curves = layer_1_curves[c(which(layer_1_curves$stroke == base_col),which(layer_1_curves$stroke != base_col)),]
  
  
  coords_L = cbind(id = path_locs_L, make_svg_df_straight(path_locs_L))
  layer_1_lines = dplyr::full_join(layer_1[layer_1$id %in% path_locs_L,], coords_L, by="id")
  layer_1_lines = layer_1_lines[c(which(layer_1_lines$stroke == base_col),which(layer_1_lines$stroke != base_col)),]
  
  
  ## Layer 2: points ##
  layer_2_points = gsub("'","",grep("ellipse", svg_xml[layer_locations$layer_loc[2]:layer_locations$layer_loc[3]], value = TRUE))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",3))))
  rx = as.numeric(gsub("'", "",gsub('"',"",(gsub('rx=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",4)))))))
  ry = as.numeric(gsub("'", "",gsub('"',"",(gsub('ry=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",5)))))))
  cx = as.numeric(gsub("'", "",gsub('"',"",(gsub('cx=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",7)))))))
  cy = as.numeric(gsub("'", "",gsub(">(.*)","",gsub('"',"",(gsub('cy=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",8))))))))
  fill = gsub("'", "",gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_2_points, "[ ]+"), "[[",6))))))
  layer_2 = data.frame(type = "ellipse", opacity, rx,ry,cx,cy,fill)
  
  ## Layer 3: Path text ##
  layer_3_text = gsub("'","",grep("text", svg_xml[layer_locations$layer_loc[3]:layer_locations$layer_loc[4]], value = TRUE))
  font_size = as.numeric(gsub('px;opacity:',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",3))))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",4))))
  fill = gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",5)))))
  x = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('x=',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",6)))))))
  y = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('y=',"",unlist(lapply(stringr::str_split(layer_3_text, "[ ]+"), "[[",7)))))))
  text = gsub("[ ]$","",gsub("</text", "",unlist(lapply(stringr::str_split(layer_3_text, ">"), "[[", 2))))
  layer_3 = data.frame(type="text", font_size, opacity, fill, x,y,text)
  
  ## layer 4: textrects ##
  layer_4_rect = gsub("'","",grep("rect", svg_xml[layer_locations$layer_loc[4]:layer_locations$layer_loc[5]], value = TRUE))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",3))))
  fill = gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",4)))))
  stroke_width = as.numeric(gsub('"',"",(gsub('stroke-width=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",5))))))
  x = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('x=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",6)))))))
  y = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('y=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",7)))))))
  height = as.numeric(gsub('"',"",(gsub('height=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",8))))))
  width = as.numeric(gsub('"',"",(gsub('width=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",9))))))
  rx = as.numeric(gsub('"',"",(gsub('rx=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",10))))))
  ry = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('ry=',"",unlist(lapply(stringr::str_split(layer_4_rect, "[ ]+"), "[[",11)))))))
  layer_4 = data.frame(type="rect", opacity, fill, stroke_width, x, y, height, width, rx,ry)
  
  ## Layer 5: Big text ##
  layer_5_text = gsub("'","",grep("text", svg_xml[layer_locations$layer_loc[5]:length(svg_xml)], value = TRUE))
  layer_5_text = layer_5_text[which(stringr::str_sub(layer_5_text, 2,5) == "text")]
  font_size = as.numeric(gsub('px;opacity:',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",3))))
  opacity = as.numeric(gsub(';',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",4))))
  fill = gsub('"',"",(gsub('fill=',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",5)))))
  x = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('x=',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",6)))))))
  y = as.numeric(gsub(">(.*)","",gsub('"',"",(gsub('y=',"",unlist(lapply(stringr::str_split(layer_5_text, "[ ]+"), "[[",7)))))))
  text = gsub("[ ]$","",gsub("<(.*)","", gsub(">(.*)", "",unlist(lapply(stringr::str_split(layer_5_text, ">"), "[[", 2)))))
  layer_5 = data.frame(type="text", font_size, opacity, fill, x,y,text)
  
  all_layers = list()
  all_layers[['layer_1_curves']] = layer_1_curves
  all_layers[['layer_1_lines']] = layer_1_lines
  all_layers[['layer_2']] = layer_2
  all_layers[['layer_3']] = layer_3
  all_layers[['layer_4']] = layer_4
  all_layers[['layer_5']] = layer_5
  
  return(all_layers)
  
}

#' Convert a svg path specification to data.frame compatible with ggplot
#' Used for paths with curves (C)
#' @param input string of svg path
#' @param rebase set the minimum path value to 0?
#' @export
#' @import stringr
#' @importFrom plyr arrange
#' @importFrom reshape2 melt
#' @importFrom plyr rbind.fill
#' @author Beth Signal
#' @examples
#'
make_svg_df = function(input, rebase = F){
  # SVG path vector types
  types = unlist(str_split(gsub("[0-9,. ]*","",input), ""))
  # locations of each path type indicator in input
  locs = unlist(lapply(str_locate_all(input, unique(types)), function(x) x[,1]))
  # make data.frame of each vector type in input, and its location within the input string
  types = str_sub(input, locs, locs)
  type_locs = rbind(data.frame(types, locs), data.frame(types = "end",locs = nchar(input)))
  type_locs = plyr::arrange(type_locs, locs)
  type_locs = type_locs[!duplicated(with(type_locs, paste(types,locs))),]
  
  # get SVG coordniates for each path type
  all_coords=NULL
  for(i in 1:(nrow(type_locs)-1)){
    all_coords = rbind(all_coords,
                       cbind(type = type_locs$types[i],
                             as.data.frame(matrix(as.numeric(unlist(str_split(
                               gsub(" ","",str_sub(input, type_locs$locs[i]+1, type_locs$locs[i+1]-1)),
                               ","))), ncol=6, byrow = TRUE))
                       )
    )
  }
  
  all_coords$type = as.character(all_coords$type)
  all_coords[which(all_coords$type != "C"),c(4:7)] = NA
  
  if(all_coords$type[nrow(all_coords)] == "Z"){
    all_coords[nrow(all_coords),-1] = all_coords[1,-1]
  }
  
  if(any(all_coords$type == "C")){
    all_coords[all_coords$type=="C",] = all_coords[all_coords$type=="C",c(1,6,7,2,3,4,5)]
  }
  
  # if TRUE start coordinates at 0
  if(rebase == TRUE){
    min_x = min(all_coords$V1)
    min_y = min(all_coords$V2)
    all_coords = all_coords
    all_coords$V1 = all_coords$V1 - min_x
    all_coords$V2 = all_coords$V2 - min_y
    all_coords$V3 = all_coords$V3 - min_x
    all_coords$V4 = all_coords$V4 - min_y
    all_coords$V5 = all_coords$V5 - min_x
    all_coords$V6 = all_coords$V6 - min_y
  }
  
  # convert from x-y points to x1,y1,x2,y2 segments (or curves)
  all_coords_plot = all_coords[-1,]
  colnames(all_coords_plot)[-1] = c("x2","y2", "xb1","yb1","xb2","yb2")
  all_coords_plot = cbind(x1=all_coords$V1[-nrow(all_coords)],y1=all_coords$V2[-nrow(all_coords)],all_coords_plot)
  all_coords_plot = all_coords_plot[,c(3,1,2,4:9)]
  # fix for Z
  all_coords_plot$type = ifelse((all_coords_plot$type == "C"), "curve","segment")
  
  if(any(all_coords_plot$type == "curve")){
    all_coords_plot_circles = all_coords_plot[all_coords_plot$type == "curve",]
    all_coords_plot_circles$bez_id = paste0(input,".b.curve", c(1:nrow(all_coords_plot_circles)))
    all_coords_plot_bez = cbind(reshape2::melt(all_coords_plot_circles[,c(2,4,6,8,10)], id.vars = c("bez_id")),
                                reshape2::melt(all_coords_plot_circles[,c(3,5,7,9,10)], id.vars = c("bez_id"))[,3])
    colnames(all_coords_plot_bez) = c("bez_id", "b_type", "x1","y1")
    all_coords_plot_bez$b_type = factor(all_coords_plot_bez$b_type, levels =c("x1","xb1","xb2","x2"))
    all_coords_plot_bez = plyr::arrange(all_coords_plot_bez, bez_id, b_type)
    all_coords_plot_bez$type="curve"
    all_coords_plot = all_coords_plot[all_coords_plot$type == "segment", c(1:5)]
    
    all_coords_plot = plyr::rbind.fill(all_coords_plot, all_coords_plot_bez)
  }else{
    all_coords_plot = all_coords_plot[, c(1:5)]
    all_coords_plot$bez_id = NA
    all_coords_plot$b_type = NA
  }
  all_coords_plot$id = input
  
  return(all_coords_plot)
}

#' Convert a svg stright line path specification to data.frame compatible with ggplot
#' Used for paths without curves (C)
#' @param input string of svg path
#' @export
#' @import stringr
#' @import ggplot2
#' @import ggforce
#' @author Beth Signal
#' @examples
#'
make_svg_df_straight = function(input){
  
  type = "segment"
  x1 = as.numeric(str_sub(unlist(lapply(str_split(input, ","),"[[",1)),2,-1))
  y1 = as.numeric(unlist(lapply(str_split(lapply(str_split(input, ","),"[[",2), "L"),"[[",1)))
  x2 = as.numeric(unlist(lapply(str_split(lapply(str_split(input, ","),"[[",2), "L"),"[[",2)))
  y2 = as.numeric(unlist(lapply(str_split(input, ","),"[[",3)))
  
  all_coords_plot = data.frame(type,x1,y1,x2,y2,curvature=NA)
  return(all_coords_plot)
}

#' make a ipath ggplot
#' @param plot_layers list of svg layer data produced by svg_to_table()
#' @export
#' @import stringr
#' @import ggplot2
#' @import ggforce
#' @author Beth Signal
#' @examples
#'
make_ipath_ggplot = function(plot_layers){
  
  layer_1_curves = plot_layers[['layer_1_curves']]
  layer_1_lines = plot_layers[['layer_1_lines']]
  layer_2 = plot_layers[['layer_2']]
  layer_3 = plot_layers[['layer_3']]
  layer_4 = plot_layers[['layer_4']]
  layer_5 = plot_layers[['layer_5']]
  
  bg_stroke = table(c(as.character(layer_1_curves$stroke), as.character(layer_1_lines$stroke), as.character(layer_2$fill)))
  bg_stroke_col = names(which.max(bg_stroke))
  
  
  layer_1 = plyr::rbind.fill(layer_1_curves, layer_1_lines)
  stroke_cols = as.character(sort(unique(layer_1$stroke)))
  layer_1_bottom = layer_1[layer_1$stroke == bg_stroke_col,]
  layer_1_top = layer_1[layer_1$stroke != bg_stroke_col,]
  
  layer_1_bottom$stroke = factor(layer_1_bottom$stroke, levels = stroke_cols)
  layer_1_top$stroke = factor(layer_1_top$stroke, levels = stroke_cols)
  
  fill_cols = as.character(sort(unique(layer_4$fill)))
  layer_4$fill = factor(layer_4$fill, levels=fill_cols)
  
  rects_left = layer_4
  rects_left$x = rects_left$x
  rects_left$width = rects_left$height
  
  rects_right = layer_4
  rects_right$x = rects_right$x + rects_right$width - (rects_right$rx*2)
  rects_right$width = rects_right$height
  rects_rounded = rbind(rects_left, rects_right)
  
  p =
    ggplot() +
    geom_bezier(aes(x = x1, y = y1, group=bez_id,col=stroke, size=stroke_width),
                data = layer_1_bottom[layer_1_bottom$type == "curve",]) +
    geom_segment(data=layer_1_bottom[layer_1_bottom$type == "segment", ],
                 aes(x=x1, y=y1, xend = x2, yend=y2, col=stroke, size=stroke_width)) +
    geom_bezier(aes(x = x1, y = y1, group=bez_id,col=stroke, size=stroke_width),
                data = layer_1_top[layer_1_top$type == "curve",]) +
    geom_segment(data=layer_1_top[layer_1_top$type == "segment", ],
                 aes(x=x1, y=y1, xend = x2, yend=y2, col=stroke, size=stroke_width)) +
    
    geom_point(data=layer_2, aes(x=cx,y=cy), size= 0.5, col=layer_2$fill[1]) +
    geom_text(data=layer_3, aes(x=x,y=y, label = text), size = 1.6, hjust = 0, vjust=0) +
    #geom_rect(data=layer_4, aes(xmin=x, ymin = y, xmax =x+width, ymax=y+height, fill=fill)) +
    geom_ellipse(data = rects_rounded, aes(x0=x+(width/2),y0=y+(height/2), a=width/2, b=height/2, angle=0, m1=2, m2=2, fill=fill), col=NA) +
    geom_rect(data = layer_4, aes(xmin=x+rx, ymin=y,ymax=y+height, xmax=x+width-rx, fill=fill), col=NA) +
    geom_text(data=layer_5, aes(x=x,y=y, label = text), size = 2.2, hjust = 0, vjust=0, fontface = "bold", col="white") +
    scale_y_continuous(trans="reverse") +
    theme_void() +
    theme(legend.position = "none") +
    scale_color_manual(values = stroke_cols)+
    scale_fill_manual(values = fill_cols) +
    scale_size_continuous(range = c(0.5,2), limits = c(2,10))
  
  return(p)
  
}

#' Change color of boxes surrounding major pathway names
#' @param pathway_map text containing svg source code
#' @param box_color hex value for box color
#' @param text_color hex value for text (in box) color
#' @return pathway map with altered box colors
#' @keywords internal
#' @import stringr
#' @author Beth Signal
recolor_textbox = function(pw_map, box_color = "#969696", text_color = '#FFFFFF'){
  
  # recolor rectangles around text
  if(!is.null(box_color)){
    g_splits = unlist(stringr::str_split(pw_map, "<g>"))
    rects = which(stringr::str_sub(g_splits,3,6) == "rect")
    g_splits[rects] = gsub("#[A-G,0-9]{6}", box_color, g_splits[rects])
  }
  #recolor text
  if(!is.null(text_color)){
    layer_locs = find_layer_locs(g_splits)
    rect_text = g_splits[layer_locs$indices[5]:layer_locs$indices[6]]
    g_splits[rect_text] = gsub("#[A-G,0-9]{6}", text_color, g_splits[rects])
  }
  pw_map2 = paste0(g_splits, collapse = "<g>")
  return(pw_map2)
}

#' Get indices for the start of each layer in svg text
#' @param g_splits vector of svg text split by the delimiter '<g>'
#' @return data.frame of layer numbers and their starting index
#' @keywords internal
#' @import stringr
#' @author Beth Signal
find_layer_locs = function(g_splits){
  
  indices = grep("layer", g_splits)
  layer_text = g_splits[indices]
  layer_number = unlist(lapply(stringr::str_split(lapply(stringr::str_split(layer_text, "layer"),"[[", 2), "'"), "[[",1))
  
  layer_number = c(layer_number, "end")
  indices = c(indices, length(g_splits))
  
  return(data.frame(layer_number, indices))
  
}

#' Reorder paths so selected (colored) paths are on top
#' @param svg_text text containing svg source code
#' @param ipath_data data.frame formatted for sending to ipath. Needs a 'color' column
#' @return pathway map (svg formatted text) with colored paths on top
#' @keywords internal
#' @import stringr
#' @author Beth Signal
colour_on_top = function(svg_text, ipath_data){
  g_splits = unlist(stringr::str_split(svg_text, "<g>"))
  layer_locs = find_layer_locs(g_splits)
  mapped_cols = unique(ipath_data$color)
  
  g_splits_by_layer = list()
  g_splits_layer = g_splits[layer_locs$indices[-nrow(layer_locs)]]
  for(j in 1:(nrow(layer_locs)-1)){
    
    if(j == (nrow(layer_locs)-1)){
      g_splits_by_layer[[j]] = g_splits[(layer_locs$indices[j]+1):(layer_locs$indices[j+1])]
    }else{
      g_splits_by_layer[[j]] = g_splits[(layer_locs$indices[j]+1):(layer_locs$indices[j+1]-1)]
    }
  }
  
  for(j in 1:(nrow(layer_locs)-1)){
    
    indices = unlist(lapply(mapped_cols, function(x) grep(x, g_splits_by_layer[[j]])))
    
    if(length(indices) > 0){
      g_split_move = g_splits_by_layer[[j]][indices]
      g_split_rm = g_splits_by_layer[[j]][-indices]
      g_split_ins = c(g_split_rm, g_split_move)
      
      g_splits_by_layer[[j]] = c(g_split_rm, g_split_move)
    }
    
    
  }
  
  g_splits_v2 = vector()
  for(j in 1:(nrow(layer_locs)-1)){
    g_splits_v2 = c(g_splits_v2, g_splits_layer[j], g_splits_by_layer[[j]])
  }
  
  svg_text2 = paste0(g_splits_v2, collapse = "<g>")
  return(svg_text2)
  
}
#' Format input data to plot with ipath
#' @param data input data.frame
#' @param id_column name of the column containing ids to plot on the ipath map. Must be one of the supported data types and in the correct format (see: https://pathways.embl.de/help.cgi)
#' @param color_column name of the column containing values to vary colors in the ipath plot
#' @param width_column name of the column containing values to vary width in the ipath plot
#' @param opacity_column name of the column containing values to vary opacity in the ipath plot
#' @param min_path_width minimum width of paths with a match
#' @param max_path_width maximum width of paths with a match. Default width for when no wwidth_column is provided.
#' @param color_type type of color palette to map values to (discrete or continous)
#' @param color_cutoff value which seperates values into higher or lower for discrete color mapping
#' @param palette RColorBrewer pallete name
#' @return data.frame formatted for sending to ipath
#' @export
#' @author Beth Signal
create_ipath_data = function(data,
                             id_column,
                             color_column=NA,
                             width_column=NA,
                             opacity_column=NA,
                             min_path_width = 3,
                             max_path_width = 10,
                             color_type = "discrete",
                             color_cutoff = 0,
                             pallete = "RdBu"){
  data = as.data.frame(data)
  output_data = as.data.frame(data[,which(colnames(data) == id_column)])
  
  if(is.character(color_column)){
    output_data = cbind(output_data,
                        color = map_values_to_colors(as.numeric(data[,which(colnames(data) == color_column)]),
                                                     color_type = color_type,color_cutoff=color_cutoff,pallete = pallete))
  }else{
    output_data = cbind(output_data,
                        color = map_values_to_colors(rep(1, nrow(output_data)),
                                                     color_type = 'discrete',color_cutoff=color_cutoff,pallete = pallete))
  }
  if(!is.na(width_column)){
    output_data$width = map_values_to_range(as.numeric(data[,which(colnames(data) == width_column)]),
                                            min_val = 5, max_val = 10)
    output_data$width = paste0("W", round(output_data$width, 2))
  }else{
    output_data$width = paste0("W", min_path_width)
  }
  if(!is.na(opacity_column)){
    output_data = cbind(output_data,
                        opacity = map_values_to_range(data[,which(colnames(data) == opacity_column)],
                                                      min_val = 0, max_val = 1))
  }else{
    output_data = cbind(output_data, opacity = 1)
  }
  
  colnames(output_data)[1] = id_column
  return(output_data)
}

#' Map a vector of values to hex colors
#' @param values numeric vector of values
#' @param color_type type of color palette to map values to (discrete or continous)
#' @param color_cutoff value which seperates values into higher or lower for discrete color mapping
#' @param palette RColorBrewer pallete name, viridis palette name (if continous), or two colors to create a colorRampPalette from.
#' @return vector of hex colors
#' @keywords internal
#' @importFrom RColorBrewer brewer.pal
#' @author Beth Signal
#' @examples
#' map_values_to_colors(seq(-4, 5))
#' map_values_to_colors(seq(-4, 5), color_type = "continous")
#' map_values_to_colors(seq(-4, 5), color_type = "continous", pallete = "magma")
#' map_values_to_colors(seq(-4, 5), color_type = "continous", pallete = c("lightblue", "blue"))

map_values_to_colors = function(values,
                                color_type = "discrete",
                                color_cutoff = 0,
                                pallete = "RdBu"){
  
  if(color_type == "discrete"){
    # make RColorBrewer pallete from pallete name
    # uses n=5 to get good insense extremes
    colors = RColorBrewer::brewer.pal(5, pallete)[c(1,5)]
    value_cols = ifelse(values > color_cutoff, colors[1], colors[2])
    return(value_cols)
  }else{
    
    ## Use n equally spaced breaks to assign each value to n-1 equal sized bins
    ii = cut(values, breaks = seq(min(values), max(values), len = 100),
             include.lowest = TRUE)
    if(pallete[1] %in% c("A", "B","C","D","E", "magma","inferno", "plasma","viridis", "cividis")){
      value_cols = viridis(100, option = pallete[1])
      
      if(all(nchar(value_cols) == 9)){
        value_cols = gsub("FF$", "", value_cols)
      }
    }else{
      if(length(pallete) == 2){
        value_cols = colorRampPalette(pallete)(99)
      }else{
        colors = RColorBrewer::brewer.pal(5, pallete)
        value_cols = colorRampPalette(colors)(99)
      }
    }
    
    ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
    value_cols = value_cols[ii]
    
  }
  
}

#' Map a vector of values to a prespecified range
#' @param values numeric vector of values
#' @param min_val value to map minimum input value to
#' @param max_val value to map maximum input value to
#' @keywords internal
#' @return vector of values
#' @author Beth Signal
#' @examples
#' map_values_to_range(seq(0, 10), min_val = 2, max_val=4)
#'
map_values_to_range = function(values, min_val= 3, max_val=10){
  
  new_values = (values - min(values))/(max(values) - min(values))
  
  new_values = (new_values*(max_val - min_val)) + min_val
  
  return(new_values)
}

#' Format ipath data.frame as parameters
#' @param ipath_data data.frame formatted for sending to ipath
#' @param default_color Hex color value. Non-selected parts of the map will use this color by default.
#' @param default_width Default path/edge width (px) for the non-selected parts of the map.
#' @param default_radius Compound radius (px) for the non-selected parts of the map.
#' @param default_opacity Opacity setting for the non-selected parts of the map. Zero represents full transparency, while one is fully opaque.
#' @param keep_colors Keep original colors? When TRUE, non-selected parts of the map will not change color to the default one.
#' @param whole_modules Select whole modules? When TRUE, any KEGG module with at least one matching edge or compound will be highlighted.
#' @param whole_pathways Select whole pathways? If TRUE, any pathway with at least one matching edge or compound will be highlighted.
#' @param query_reactions Query reaction compounds. If TRUE, compound presence within each edges reactions will also be checked
#' @param tax_filter string containing either NCBI tax ID or KEGG 3 letter species code(s). Only pathways present in selected species will be included in the map.
#' @param export_dpi numeric value for DPI of SVG output.
#' @keywords internal
#' @return list of parameters for ipath
#' @author Beth Signal
to_parameters = function(ipath_data,
                         default_color = '#bdbdbd',
                         default_width = 3,
                         default_radius = 7,
                         default_opacity = 1,
                         keep_colors = FALSE,
                         whole_modules = FALSE,
                         whole_pathways = FALSE,
                         query_reactions = FALSE,
                         tax_filter = '',
                         export_dpi = 1200){
  
  selection  = paste0(apply(ipath_data, 1, function(x) paste(x, collapse =" ")),
                      collapse = "\n")
  
  ipath_parameters = list(selection = selection,
                          export_type = 'svg',
                          keep_colors = ifelse(keep_colors, 1, 0),
                          include_metabolic = 1,
                          include_secondary = 0,
                          include_antibiotic = 0,
                          include_microbial = 0,
                          whole_modules = ifelse(whole_modules, 1, 0),
                          default_opacity = ifelse(default_opacity, 1, 0),
                          whole_pathways = ifelse(whole_pathways, 1, 0),
                          default_width = default_width,
                          default_color = default_color,
                          default_radius = default_radius,
                          query_reactions = ifelse(query_reactions, 1, 0),
                          tax_filter = tax_filter,
                          export_dpi = export_dpi
  )
  return(ipath_parameters)
}

try_match = function(input_a, input_b, output ="inputs_short"){
  
  
  if(length(grep(input_a[1], input_b)) > 0){
    # input b is longer
    longer = "b"
    input_long = input_b
    input_short = input_a
    
  }else if(length(grep(input_b[1], input_a)) > 0){
    # input a is longer
    longer = "a"
    input_long = input_a
    input_short = input_b
  }else{
    stop(paste0("Cannot find ", input_a[1], " or ", input_b[1], " within paired list."))
  }
  
  # input a is longer
  index = grep(input_short[1], input_long)
  
  
  location = str_locate(input_long[index][1],input_short[1])
  end = location[1,2]
  start = location[1,1]
  if(end < nchar(input_long[index][1])){
    
    ncount = 2
    plus_n = 1
    while(ncount > 1 & end+plus_n <= nchar(input_long[index][1])){
      subseq = str_sub(input_long[index][1], end+1, end+plus_n)
      ncount = str_count(input_long[index][1], subseq)
      plus_n = plus_n + 1
    }
    input_long = unlist(lapply(str_split(input_long, subseq),"[[",1))
  }
  
  if(start > 1){
    
    ncount = 2
    plus_n = 0
    while(ncount > 1 & start+plus_n <= nchar(input_long[index][1])){
      subseq = str_sub(input_long[index][1], start, start+plus_n)
      ncount = str_count(input_long[index][1], subseq)
      plus_n = plus_n + 1
    }
    input_long = unlist(lapply(str_split(input_long, subseq),"[[",2))
  }
  
  if(longer == "a"){
    m = match(input_long, input_short)
    input_a = input_long
    input_b = input_short
  }else{
    m = match(input_short, input_long)
    input_a = input_short
    input_b = input_long
  }
  if(output == "index"){
    return(m)
  }else{
    new_values = list()
    new_values[[1]] = input_a
    new_values[[2]] = input_b
    return(new_values)
  }
}

join_and_drop = function(data1, data2, col1, col2){
  
  data1$join_id = data1[,match(col1, colnames(data1))]
  data2$join_id = data2[,match(col2, colnames(data2))]
  
  data1[,match(col1, colnames(data1))] = NULL
  data2[,match(col2, colnames(data2))] = NULL
  
  data = dplyr::left_join(data1, data2, by = 'join_id')
  data = cbind(first_col = data$join_id, data)
  data = data[,-which(colnames(data) == 'join_id')]
  colnames(data)[1] = col1
  
  return(data)
}

# Define UI for application that draws a histogram
ui = navbarPage("Metabolic Pathway Visualisation",
                
                tabPanel("Intro",
                         includeMarkdown("./md/intro.md"),
                         hr()),
                
                tabPanel("Upload data",
                         sidebarLayout(
                           sidebarPanel(
                             # upload files
                             fileInput("file1", "Upload CSV file with Pathway IDs",
                                       multiple=FALSE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             radioButtons("file1head", "Use first row as column names?",
                                          choices=c("Yes", "No"), selected = "Yes"),
                             fileInput("file2", "[Optional] Upload CSV file of corresponding gene/protein differential expression values",
                                       multiple=FALSE,
                                       accept = c("text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")),
                             radioButtons("file2head", "Use first row as column names?",
                                          choices=c("Yes", "No"), selected = "Yes"),
                             # select geneids for each of the two datasets
                             uiOutput("var_ui_geneid1"),
                             uiOutput("var_ui_geneid2")
                             
                           ),
                           mainPanel(
                             tabPanel('Input data',tableOutput('input_data'))
                           )
                         )
                ),
                
                tabPanel("View pathways",
                         sidebarLayout(
                           sidebarPanel(
                             uiOutput("var_ui_id"),
                             uiOutput("var_ui_col"),
                             uiOutput("var_ui_width"),
                             
                             selectInput("color_type","Colour type:", choices = c("continous", "discrete"), selected = "discrete"),
                             selectInput("color_pal","Colour palette:", choices = c("magma","inferno", "plasma","viridis", "cividis"), selected = "magma"),
                             #radioButtons("showtext2", "Show pathway names?", choices = c("Yes", "No"), selected="Yes"),
                             #radioButtons("showtext1", "Show pathway group names?", choices = c("Yes", "No"), selected="Yes"),
                             
                             downloadButton("download_svg", "Download SVG"),
                             downloadButton("download_pdf", "Download PDF")
                           ),
                           mainPanel(
                             plotOutput("pathPlot")
                           )
                         )
                         
                )
                
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  data <- reactive({
    req(input$file1)
    head = ifelse(input$file1head=="Yes", T, F)
    message(tools::file_ext(input$file1$name))
    message(input$file1$name)
    if(tools::file_ext(input$file1$name) == "csv"){
      data = read.csv(input$file1$datapath, header=head)
    }else{
      data = read.delim(input$file1$datapath, header=head)
    }
    
    data_p = data
    message("Uploaded")
    
    if(!is.null(input$file2)){
      message("Passed")
      message(str(input$file2))
      
      functional_annot=data
      
      expression_metadata = data_e()
      col1 = match(input$geneid1_col, colnames(functional_annot))
      col2 = match(input$geneid2_col, colnames(expression_metadata))
      m = try_match(functional_annot[,col1], expression_metadata[,col2])
      functional_annot[,col1] = m[[1]]
      expression_metadata[,col2] = m[[2]]
      
      data_p = join_and_drop(functional_annot, expression_metadata, input$geneid1_col, input$geneid2_col)
      
    }else{
      message("Skipped")
      data_p = data
    }
    
    data_p
  })
  
  data_e <- reactive({
    req(input$file2)
    req(input$file1)
    if(tools::file_ext(input$file2) == "csv"){
      data_e = read.csv(input$file2$datapath)
    }else{
      data_e = read.delim(input$file2$datapath)
    }
    data_e
    
  })
  
  observe({
    
    col_types = input$color_type
    if(col_types == "discrete"){
      updateSelectInput(session, "color_pal","Colour palette:",
                        choices = c("RdBu", "PuOr","PRGn","PiYG","BrBG"), selected = "RdBu")
    }else{
      updateSelectInput(session, "color_pal","Colour palette:",
                        choices = c("magma","inferno", "plasma","viridis", "cividis"), selected = "magma")
    }
    
  })
  
  # data_p <- reactive({
  #   req(data())
  #   #req(input$geneid1_col)
  #   #req(input$geneid2_col)
  #
  #   functional_annot = data()
  #
  #
  #   if(!is.null(data_e()) & !is.null(input$geneid2_col)& !is.null(input$geneid1_col)){
  #
  #     expression_metadata = data_e()
  #     col1 = match(input$geneid1_col, colnames(functional_annot))
  #     col2 = match(input$geneid2_col, colnames(expression_metadata))
  #     m = try_match(functional_annot[,col1], expression_metadata[,col2])
  #     functional_annot[,col1] = m[[1]]
  #     expression_metadata[,col2] = m[[2]]
  #
  #     data_p = join_and_drop(functional_annot, expression_metadata, input$geneid1_col, input$geneid2_col)
  #
  #   }else{
  #     data_p = data()
  #   }
  #
  #   data_p
  # })
  
  
  
  svg_xml <- reactive({
    req(input$file1)
    req(input$color_col)
    req(input$width_col)
    
    if(input$color_col == "no color"){
      color_col=NA
    }else{
      color_col = input$color_col
    }
    
    message(color_col)
    
    if(input$width_col == "no width"){
      width_col=NA
    }else{
      width_col = input$width_col
    }
    
    if(input$id_col %in% colnames(data())){
      
      ipath_data = create_ipath_data(data(),
                                     id_column = input$id_col,
                                     color_column = color_col,
                                     color_type = input$color_type,
                                     pallete = input$color_pal,
                                     width_column = width_col)
      message("create_ipath_data")
      
      svg_xml = print_ipath_pdf(ipath_data = ipath_data)
      message("print_ipath_data")
    }
    svg_xml
  })
  
  
  output$var_ui_id <- renderUI({
    id_options = colnames(data())
    
    default_id_col = ifelse("id" %in% id_options | length(id_options) == 0,  "id", id_options[2])
    id_options = id_options[which(id_options != "id")]
    
    
    
    selectInput("id_col", "Choose variable containing pathway IDs:", choices= c("id", id_options), selected = default_id_col)
  })
  
  output$var_ui_col <- renderUI({
    selectInput("color_col", "Choose variable to color pathways:", choices= c("no color", colnames(data())), selected = "no color")
  })
  
  output$var_ui_width <- renderUI({
    selectInput("width_col", "Choose variable to change width of pathways:", choices= c("no width", colnames(data())), selected = "no width")
  })
  
  output$var_ui_geneid1 <- renderUI({
    selectInput("geneid1_col", "Choose column name with sequence ids in first dataset:", choices= c(colnames(data())), selected = colnames(data())[1])
  })
  output$var_ui_geneid2 <- renderUI({
    selectInput("geneid2_col", "Choose column name with sequence ids in second dataset:", choices= c("NA", colnames(data_e())), selected = colnames(data_e())[1])
  })
  
  output$download_svg <- downloadHandler(
    
    filename = function() {
      gsub("[.]txt", ".svg", gsub("[.]csv", ".svg", input$file1))
    },
    content = function(file) {
      write.table(svg_xml(), file, row.names = FALSE, col.names = FALSE, quote=F, sep="\n")
    },
    contentType = "svg"
    
  )
  
  output$download_pdf <- downloadHandler(
    
    filename = function() {
      gsub("[.]txt", ".pdf", gsub("[.]csv", ".pdf", input$file1))
    },
    content = function(file) {
      ggsave(file, pathway_plot(), device = 'pdf',width=14,height=8.5)
    },
    contentType = "pdf"
    
  )
  
  pathway_plot = reactive({
    plot_layers = svg_to_table(svg_xml())
    p = make_ipath_ggplot(plot_layers = plot_layers)
    p
  })
  
  output$pathPlot <- renderPlot({
    
    print(pathway_plot())
    
  })
  
  output$input_data <- renderTable({
    data()
  })
  
  # output$contents <- renderTable({
  #
  #     req(input$file1)
  #     uniprot_combined <- read.csv(input$file1$datapath)
  #     ipath_data = create_ipath_data(uniprot_combined,
  #                                    id_column = input$id_column)
  #     svg_xml = print_ipath_pdf(ipath_data = ipath_data)
  #     plot_layers = svg_to_table(svg_xml)
  #     p = make_ipath_ggplot(plot_layers = plot_layers)
  #
  #     return(head(uniprot_combined))
  # })
}


# Run the application
shinyApp(ui = ui, server = server)
