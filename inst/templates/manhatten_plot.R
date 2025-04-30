# Load required libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(ggrepel)
library(ggtext)
# source('scripts/util.R')
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))
# specify colour scheme
SEQUOIA_COLS = c("#00A79D","#E67E25","#F0C418", "#bc5d92","#885dbc","#c01d41",
                 "#a4d3db","#62BC5D","#074da3","#8DC63F", "#287D7D","#1C344D",
                 "#23B99A","#69dafe","#e570ac","#4dbd8e")


# load data
df = read_excel('results/id_exposure_rs12199003_unscaled_Summary.xlsx',sheet='Clinical_FinnGen_r10')

# create a direction column
df$dir=as.factor(
  ifelse(
    df$Pfdr<=0.05 & df$beta<0,
    "neg_sig",
    ifelse(
      df$Pfdr>0.05 & df$beta<0, 
      "neg_nonsig",
      ifelse(
        df$Pfdr<=0.05 & df$beta>0, 
        "pos_sig", 
        "pos_nonsig"
      )
    )
  )
)

# create label_group column
df$label_group = factor(
  df$label_group, levels = c(sort(setdiff(unique(df$label_group), "Other")), "Other")
)

# Order the data frame based on label_group and trait
df = df[order(df$label_group, df$trait, decreasing = F),]

# create helper columns
df$sig=as.factor(ifelse(df$Pfdr<=0.05, TRUE, FALSE))

# specify colour for plotting
df$coluse=SEQUOIA_COLS[df$label_group]

# create label group label with colour
df$plotlabel=paste("<span style = 'color: ",
                   df$coluse,";'>",df$label_group,
                   "</span>", sep = "")
df$plotlabel=as.factor(df$plotlabel)


# create meanid for label position
mygroup = c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150)
df$groupnum=mygroup[df$label_group]
df$count=1:nrow(df)
df$count2=df$count+df$groupnum
df=df %>% 
  group_by(label_group) %>%
  mutate(id=row_number()+groupnum) %>%
  mutate(maxcount=max(count2)) %>%
  mutate(meanid=mean(count2)) %>%
  ungroup()

toname=df$trait[df$Pfdr<=0.05]
ylim = round(max(-log10(df$Pfdr))+ 1) 
ylim = max(ylim, 1.4)

f1 = ggplot(df, aes(x=count2, y=-log10(Pfdr))) +
  # Show all points
  geom_point(aes(shape = as.factor(dir), color = coluse, fill=coluse, alpha = sig), size = 2) +
  scale_x_continuous(limits = c(0,max(df$count2)),breaks= df$meanid,  label = df$plotlabel, expand=c(0.02,0)) +
  scale_colour_identity()+
  scale_fill_identity()+
  scale_shape_manual(values = c("pos_sig"=24, "pos_nonsig"=24, "neg_sig"=25,"neg_nonsig"=25)) +
  geom_hline(aes(yintercept=1.30103), linetype =2, colour="black")+
  scale_y_continuous(expand = c(0, 0), limits = c(0,ylim), 
                     name = "-log<sub>10</sub>(p<sub>FDR</sub>)") + 
  geom_label_repel(data=df[df$trait%in%toname,],
                   aes(label=toname),
                   size=3, force=1.3)+
  theme_bw(base_size = 14) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    axis.text.x = element_markdown(angle = 45, hjust = 1),
    axis.title.y = element_markdown(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )+ 
  labs(caption = "Association scale per 1 SD higher genetically predicted BMI through GFRAL")

# save plot
ggsave("GDF15_phewas_manhattan_plot.png", 
       f1, dpi=300, height=180, width=330, units="mm")