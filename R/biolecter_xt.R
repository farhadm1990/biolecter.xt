biolect_xt <- function(path_to_biolector_xlsx, working_dir, caliberation = "both", replic = TRUE, n_sheets = 8){ 

setwd(working_dir)

#loading required packages
if(!requireNamespace("pacman", quietly = TRUE)){
  install.packages("pacman")
    }

    library("pacman")
    pacman::p_load(data.table, ggrepel, glue, stringr, reshape2, tidyverse, readxl, patchwork)


#Loading data and data wrangling
path_to_data = path_to_biolector_xlsx 
#importing all data sheets in a list
all_dat = list()
for(i in seq_along(1:n_sheets)){
all_dat[[i]] = readxl::read_xlsx(path = path_to_data, sheet =i)

    }
layout_df = all_dat[[2]][-c(1:3),]
#metadata
metdat = all_dat[[1]] %>% as.matrix()
user = metdat[1,2]
date = date()

#process data
process_data = all_dat[[3]]

#retrieving channel names
n_chan  = n_sheets - 4 
chans_name = paste("Channel", 1:n_chan, "Name")
chans_gain = paste("Channel", 1:n_chan, "Gain")

filterset = c()
for(chan in chans_name){
filterset[[chan]] = metdat[metdat[,1] == chan,2][[1]]
}
filterset = gsub(" ", "_", filterset)

#Extracting gain values
chans_gain_value = c()

for(gain in chans_gain){
chans_gain_value[[gain]] <- metdat[metdat[,1] == gain,2][[1]]
}


#chan_dat <- all_dat[length(chans_name):(length(all_dat)-1)]
chan_dat <- all_dat[4:(length(all_dat)-1)]

filterset_cal = paste0(filterset, "_cal")
filterset_raw = paste0(filterset, "_raw")


#Splitting dataset into clibrated and raw for each channel
##calibrated
chan_cal= list()
for(i in seq_along(1:length(chan_dat))){
conm <- colnames(chan_dat[[i]])
cyc_tm <- chan_dat[[i]][,c(1,2)]
filt_name <- filterset_cal[i]

chan_cal[[i]] <- chan_dat[[i]][, grepl(pattern = glue("Ch{i}_cal_*"), x = conm)]
final_df <- data.frame(cyc_tm, chan_cal[[i]])

colnames(final_df) <- as.character(c("cycle", 
                                    "time", 
                                    chan_cal[[i]][3,]))
final_df = final_df[-c(1:3),] %>% mutate(filterset = filt_name)

#adding filterset column

assign(x= filt_name, final_df)
}
rm(chan_cal)

#Raw
chan_raw = list()
for(i in seq_along(1:length(chan_dat))){
conm <- colnames(chan_dat[[i]])
cyc_tm <- chan_dat[[i]][,c(1,2)]
filt_name <- filterset_raw[i]

chan_raw[[i]] <- chan_dat[[i]][, grepl(pattern = glue("Ch{i}_raw_*"), x = conm)]
final_df <- data.frame(cyc_tm, chan_raw[[i]])

colnames(final_df) <- as.character(c("cycle", 
                                    "time", 
                                    chan_raw[[i]][3,]))
final_df = final_df[-c(1:3),] %>% mutate(filterset = filt_name)
assign(x= filt_name, final_df)
}

rm(chan_raw)

full = list()
for(i in 1:length(filterset)){
    cal <- get(filterset_cal[i]) #calling a dataframe by a string using get()
    raw <- get(filterset_raw[i])
    full[[i]] = rbind(cal, raw)
    
}

 
full_long = data.table::rbindlist(full) %>% 
                pivot_longer(!c("cycle", "time", "filterset"), 
                names_to = "well", 
                values_to = "value") %>% 
                filter(!is.na(value)) %>% 
                mutate(gain = ifelse(grepl(filterset, 
                pattern = "Biomass_1"), paste0(chans_gain_value[[1]]),
                ifelse(grepl(filterset, 
                pattern = "Biomass_2"), paste0(chans_gain_value[[2]]), 
                ifelse(grepl(filterset, 
                pattern = "Biomass_3"), 
                paste0(chans_gain_value[[3]]), "pH")))) %>%
                mutate(calibrated = as.factor(ifelse(grepl(filterset, pattern = "_cal"), "Calibrated", "Raw")))

full_long$value <- as.numeric(full_long$value)

full_long$filterset = str_replace_all(c("_cal"= "", "_raw" = ""), string = full_long$filterset)

#a function for creating multiplots

mutli_plot = function (data, 
                        well, 
                        flt_set = filterset, 
                        calibrated = TRUE,
                        isolate = isolate) { 

colrs = c("darkmagenta", "forestgreen", "deepskyblue", "darkorange")
#Making a folder for the user
if(dir.exists(paths =  glue("{getwd()}/{user}"))){
output_dir = glue("{getwd()}/{user}")
} else {
dir.create(glue("{getwd()}/{user}"))
output_dir = glue("{getwd()}/{user}")
}    
    
#Creating subsets of data for each fitlerset 
pattern_flt = gsub("_.*", replacement = "", flt_set)
anaerob = metdat[grepl(pattern = "anaerobic", ignore.case = TRUE, x = metdat[,1]),2]
#Biologcial replicate?

if(replic == TRUE){
    #Replicate On

    lay_df = layout_df[complete.cases(layout_df),] %>% group_by(Description) %>% data.frame() %>% mutate(group = as.numeric(factor(Description))) 
    #expanding group vairable to the length of long data
    if(length(unique(lay_df[lay_df$Well %in% data$well, "group"])) > 1){
    data$group <- rep(lay_df[lay_df$Well %in% data$well, "group"], length(data$well)/length(unique(data$well)))
    data$isolate <- rep(lay_df[lay_df$Well %in% data$well, "Description"], length(data$well)/length(unique(data$well)))
    } else if((length(unique(lay_df[lay_df$Well %in% data$well, "group"])) == 1)){
    data$group <- rep(unique(lay_df[lay_df$Well %in% data$well, "group"]), nrow(data))
    data$isolate <- rep(unique(lay_df[lay_df$Well %in% data$well, "Description"]), nrow(data))
    }
    
    #calculating mean and std
    for(i in isolate){

    df <- data %>% group_by(time, isolate, filterset, gain, calibrated) %>% summarise(mean = mean(value), std = sd(value), .groups = "drop") %>% filter(isolate == i ) 

    if(calibrated == TRUE){
    df = df %>% filter(calibrated == "Calibrated")
    calib = "calibrated"
    } else if (calibrated == FALSE){
    df = df %>% filter(calibrated == "Raw")
    calib = "raw"
    } else if (calibrated == "both") {
    df = df
    calib = "calibrated and raw"
    } else {
    stop("Please enter a valid value, TRUE, FALSE, or both for clibrated argument.")
    }

    #biomass
    df_biomass <- df %>% filter(grepl("Biomass", filterset)) %>%  mutate(gain2 = glue("Biomass (gain {gain})"), gain = NULL) 

#ph 
    df_ph <- df %>% filter(grepl("pH", filterset)) %>%  mutate( gain = NULL) 

#DO 
    df_do <- df %>% filter(grepl("DO", filterset)) %>%  mutate( gain = NULL) 

#plots
    colrsp1 = colrs[1:length(unique(df_biomass$gain2))]
    p1 = df_biomass %>% ggplot(aes(x = time, y = mean, group = gain2)) +
            geom_pointrange(aes(ymin = mean-std, ymax = mean+std), 
            color = "lightblue", alpha = 0.25 ) + 
            geom_errorbar(aes(ymin = mean-std, ymax = mean+std), 
            width = 1, color = "deeppink1", alpha = 0.05,
            position = position_dodge(0.5)) +
            geom_line( aes( color = gain2), show.legend = TRUE, lwd = 0.5) +
            geom_point(pch = 21, aes(fill = gain2), 
            color = "white", stroke = 0.05, show.legend = FALSE) +
            scale_color_manual(values = colrsp1) +
            scale_fill_manual(values = colrsp1) +
            scale_x_continuous(n.breaks = 8) + 
            theme_bw()  + 
            facet_wrap(~calibrated, scale = "free_y")  +
            labs(color = "Filterset") +
            guides(fill = FALSE) +
            geom_label_repel(data = . %>% filter(time == last(time)), 
            aes(label = gain2), color = "black", show.legend = FALSE,
            na.rm = TRUE, size = 02, 
            nudge_x = 1, 
            nudge_y = 1,
            box.padding = unit(0.005,  "line"), 
            label.size = 0.1) + 
            ylab("Biomass") +
            xlab("") + 
            theme(strip.text = element_text(face = "bold", color = "white"),
            strip.background = element_rect( fill = "aquamarine4")) + 
            ggtitle(glue("Timeserie plot of different filtersets for isolate {i} with biological replicate."), subtitle = glue("Anaerobic mode: {anaerob}"))

            #a control for empty plots
            if(dim(p1$data)[1] == 0){
            p1 <- NULL
            }
            

        colrsp2 = colrs[1:length(unique(df_ph$filterset))]
        p2 = df_ph %>% ggplot(aes(x = time, y = mean)) +
            geom_pointrange(aes(ymin = mean-std, ymax = mean+std), 
            color = "deepskyblue", alpha = 0.25 ) + 
            geom_errorbar(aes(ymin = mean-std, ymax = mean+std), 
            width = 1, color = "deeppink1", alpha = 0.05,
            position = position_dodge(0.5)) +
            geom_line( aes(group = filterset, 
            color = filterset), show.legend = FALSE, lwd = 0.5) +
            geom_point(pch = 21, aes(fill = filterset), 
            color = "white", stroke = 0.05, show.legend = FALSE) +
            scale_color_manual(values = colrsp2) +
            scale_fill_manual(values = colrsp2) +
            scale_x_continuous(n.breaks = 8) + 
            theme_bw()  + 
            facet_wrap(~calibrated, scale = "free_y")  +
            labs(color = "Filterset")  + 
            ylab("pH") +
            xlab("Time, h") + 
            theme(strip.text = element_text(face = "bold", color = "white"),
            strip.background = element_rect( fill = "aquamarine4"))+
            ggtitle(glue("Timeserie plot of different filtersets for isolate {i} with biological replicate."), subtitle = glue("Anaerobic mode: {anaerob}"))
            #a control for empty plots
            if(dim(p2$data)[1] == 0){
            p2 <- NULL
            }
            


        #DO plot
        colrsp3 = colrs[1:length(unique(df_do$filterset))]
        if (anaerob == "Off"){
        p3  = df_do %>%   
        ggplot(aes(x = time, y = mean)) +
        geom_pointrange(aes(ymin = mean-std, ymax = mean+std), 
        color = "deepskyblue", alpha = 0.25 ) + 
        geom_errorbar(aes(ymin = mean-std, ymax = mean+std), 
        width = 1, color = "deeppink1", alpha = 0.05,
        position = position_dodge(0.5))+
        geom_line( aes(y = value, group = filterset, 
        color = filterset), show.legend = FALSE, lwd = 0.5) +
        geom_point(pch = 21, aes(fill = filterset), 
        color = "grey", stroke = 0.05, show.legend = FALSE) +
        scale_color_manual(values = colrsp3) +
        scale_fill_manual(values = colrsp3) +
        scale_x_continuous(n.breaks = 8) + 
            theme_bw()  +
            facet_wrap(~calibrated, scale = "free_y") + 
            labs(color = "Filterset") +
            geom_label_repel(data = . %>% filter(time == last(time)), 
            aes(label = filterset), color = "black", show.legend = FALSE,
            na.rm = TRUE, size = 02, 
            nudge_x = 1, 
            nudge_y = 1,
            box.padding = unit(0.005,  "line"), 
            label.size = 0.1) + 
            ylab("DO") +
            xlab("Time, h") +  
            theme(strip.text = element_text(face = "bold", color = "white"),
            strip.background = element_rect( fill = "aquamarine4"))

            if(dim(p3$data)[1] == 0){
            p3 <- NULL
            }
            if(!is.null(p1) & !is.null(p2) & !is.null(p3)){
            p3$theme$strip.text <- element_blank()
            p2$theme$strip.text <- element_blank()
            p3$theme$strip.text <- element_blank()
            p2$theme$strip.text <- element_blank()
            p3$labels$title <- NULL
            p3$labels$subtitle <- NULL
            p2$labels$title <- NULL
            p2$labels$subtitle <- NULL
            p1 + p2 + p3 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            } else if(!is.null(p1) & is.null(p2) & !is.null(p3)){
            p3$theme$strip.text <- element_blank()
            p3$labels$title <- NULL
            p3$labels$subtitle <- NULL
            p1 + p3 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            } else if (is.null(p1) & is.null(p2) & !is.null(p3)){
            p3 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            }
              
      
            } else {
            if(!is.null(p1) & !is.null(p2)){
            p2$theme$strip.text <- element_blank()
            p2$labels$title <- NULL
            p2$labels$subtitle <- NULL
            p1 + p2+ patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            } else if(!is.null(p1) & is.null(p2)){
            p1 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            }

        }
        ggsave(glue("{output_dir}/{user}_Timeseries_plot_of_{calib}_data_for_isolate_{i}_with_{length(flt_set)}_filtersets (Anaerobic mode {anaerob}). Replicate_date:{date}.jpeg"), dpi = 700)

    }#end foor loop for isolate

for(i in length(1:3)){
if(length(ls()[ls()==glue("p{i}")])>0){ 

} else {

}
}


}  else {


    for(i in 1:length(well)){
    df = data %>% filter(well == well[i] )

    if(calibrated == TRUE){
    df = df %>% filter(calibrated == "Calibrated")
    calib = "calibrated"
    } else if (calibrated == FALSE){
    df = df %>% filter(calibrated == "Raw")
    calib = "raw"
    } else if (calibrated == "both") {
    df = df
    calib = "calibrated and raw"
    } else {
    stop("Please enter a valid value, TRUE, FALSE, or both for clibrated arbument.")
    }


#Biomass df
df_biomass =  df[grepl("Biomass", df$filterset), ] %>% 
                mutate(gain2 = glue("Biomass (gain {gain})"))

#pH df
df_ph = df %>% filter(filterset == "pH")

#DO df
df_do = df %>% filter(filterset == "DO")



#Biomas plot
colrsp1 = colrs[1:length(unique(df_biomass$gain2))]
p1  = df_biomass %>%   
    ggplot(aes(x = time, y = value)) +
        geom_line( aes(y = value, group = gain2, color = gain2), lwd = 0.5, show.legend = TRUE) +
        geom_point(pch = 21, aes(fill = gain2), 
        color = "grey", stroke = 0.05, show.legend = FALSE) +
        scale_color_manual(values = colrsp1) +
        scale_fill_manual(values = c("darkmagenta", "forestgreen", "deepskyblue")) +
        scale_x_continuous(n.breaks = 8) + 
            theme_bw()  +
            facet_wrap(~calibrated, scale = "free_y") + 
            labs(color = "Filterset") +
            guides(fill = FALSE) +
            geom_label_repel(data = . %>% filter(time == last(time)), 
            aes(label = gain2), color = "black", show.legend = FALSE,
            na.rm = TRUE, size = 02, 
            nudge_x = 0.1, 
            nudge_y = 0.1,
            box.padding = unit(0.005,  "line"), 
            label.size = 0.1) + 
            ylab("Biomass for there gains") +
            xlab("")+ 
            theme(strip.text = element_text(face = "bold", color = "white"),
                 strip.background = element_rect( fill = "aquamarine4")) + 
            ggtitle(glue("Timeserie plot of different filtersets for well {well[i]}: {isolate[i]}."), subtitle = glue("Anaerobic mode: {anaerob}"))
            #a control for to skip empty plots
            if(dim(p1$data)[1] == 0){
            p1 <- NULL
            }
#ph plot
colrsp2 = colrs[1:length(unique(df_ph$filterset))]
p2  = df_ph %>%   
    ggplot(aes(x = time, y = value)) +
        geom_line( aes(y = value, group = filterset, 
        color = filterset), show.legend = FALSE, lwd = 0.5) +
        geom_point(pch = 21, aes(fill = filterset), 
        color = "grey", stroke = 0.05, show.legend = FALSE) +
        scale_color_manual(values = colrsp2) +
        scale_fill_manual(values = colrsp2) +
        scale_x_continuous(n.breaks = 8) + 
            theme_bw()  +
            facet_wrap(~calibrated, scale = "free_y") + 
            labs(color = "Filterset") + 
            ylab("pH") +
            xlab("Time, h") + 
            theme(strip.text = element_text(face = "bold", color = "white"),
            strip.background = element_rect( fill = "aquamarine4")) +
            ggtitle(glue("Timeserie plot of different filtersets for well {well[i]}: {isolate[i]}."), 
            subtitle = glue("Anaerobic mode: {anaerob}"))
            
            if(dim(p2$data)[1] == 0){
            p2 <- NULL
            }
            
#DO plot
colrsp3 = colrs[1:length(unique(df_do$filterset))]
if (anaerob == "Off"){
    p3  = df_do %>%   
    ggplot(aes(x = time, y = value)) +
        geom_line( aes(y = value, group = filterset, 
        color = filterset), show.legend = FALSE, lwd = 0.5) +
        geom_point(pch = 21, aes(fill = filterset), 
        color = "grey", stroke = 0.05, show.legend = FALSE) +
        scale_color_manual(values = colrsp3) +
        scale_fill_manual(values = colrsp3) +
        scale_x_continuous(n.breaks = 8) + 
            theme_bw()  +
            facet_wrap(~NULL, scale = "free_y") + 
            labs(color = "Filterset") +
            geom_label_repel(data = . %>% filter(time == last(time)), 
            aes(label = filterset), color = "black", show.legend = FALSE,
            na.rm = TRUE, size = 02, 
            nudge_x = 1, 
            nudge_y = 1,
            box.padding = unit(0.005,  "line"), 
            label.size = 0.1) + 
            ylab("DO") +
            xlab("Time, h") +  
            theme(strip.text = element_text(face = "bold", color = "white"),
            strip.background = element_rect( fill = "aquamarine4"))

            if(dim(p3$data)[1] == 0){
            p3 <- NULL
            }
            
            if(!is.null(p1) & !is.null(p2) & !is.null(p3)){
            p3$theme$strip.text <- element_blank()
            p2$theme$strip.text <- element_blank()
            p3$labels$title <- NULL
            p3$labels$subtitle <- NULL
            p2$labels$title <- NULL
            p2$labels$subtitle <- NULL
            p1 + p2 + p3 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            } else if(!is.null(p1) & is.null(p2) & !is.null(p3)){
            p3$theme$strip.text <- element_blank()
            p3$labels$title <- NULL
            p3$labels$subtitle <- NULL
            p1 + p3 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            } else if (is.null(p1) & is.null(p2) & !is.null(p3)){
            p3 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            }
              
      
            } else {
            if(!is.null(p1) & !is.null(p2)){
            p2$theme$strip.text <- element_blank()
            p2$labels$title <- NULL
            p2$labels$subtitle <- NULL
            p1 + p2+ patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            } else if(!is.null(p1) & is.null(p2)){
            p1 + patchwork::plot_layout(ncol =1, nrow =length(unique(pattern_flt)))
            }
            }
            ggsave(glue("{output_dir}/User: {user}_Timeseries_plot_of_{calib}_data_for_well_{well[i]}_with_{length(flt_set)}_filtersets (Anaerobic mode {anaerob})_date: {date}.jpeg"), dpi = 700)
}

}

}




                   


#Generating the plots

well = unique(full_long$well)
gains = levels(full_long$gain)
isolate = pull(layout_df[layout_df$Well %in% well, "Description"])
filt_set = filterset



#running the buildin function
mutli_plot(data= full_long, 
            well = well, 
            calibrated = caliberation, 
            flt_set = filt_set,
            isolate = isolate)

output_dir = glue("{getwd()}/{user}")

write.table(full_long, glue("{output_dir}/full_long.tsv"), sep = "\t")

write.table(metdat, glue("{output_dir}/metdat.tsv"), sep = "\t")

write.table(layout_df, glue("{output_dir}/layout_df.tsv"), sep = "\t")

write.table(process_data, glue("{output_dir}/process_data.tsv"), sep = "\t")

            print(glue("Plots and tables have been put in {working_dir}"))
            print("Thanks for using this script :)")
}