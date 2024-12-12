### load the library
library(tidyverse)
library(ggplot2)

### load the data
coverage_unfilter <- read_tsv("./window_depth.bed",
                              col_names = c("Reference","Start","End","Coverage")) %>% 
  mutate(Mode = "unfiltered")
coverage_mq20 <- read_tsv("./window_depth_mq20.bed",
                              col_names = c("Reference","Start","End","Coverage")) %>% 
  mutate(Mode = "MQ20")

coverage <- bind_rows(coverage_mq20,coverage_unfilter)

coverage %>% 
  mutate(Start = Start/10000) %>% 
  ggplot(., aes(x = Start, y = Coverage, color = Mode, fill = Mode)) +
  geom_line(size = 1) + # Draw lines
  geom_ribbon(aes(ymin = 0, ymax = Coverage), alpha = 0.3) + # Fill area under the line
  scale_x_continuous(name = "Genome Position (Kb)") +
  scale_y_continuous(name = "Coverage") +
  scale_color_manual(values=c('#999999','#E69F00'))+
  scale_fill_manual(values=c('#999999','#E69F00'))+
  theme_bw() +
  labs(color = "Mode",
       fill = "Mode") +
  theme(text = element_text(size = 12),
        legend.position = "bottom")

ggsave("./target_species.coverage.pdf",width = 10, height = 7)
