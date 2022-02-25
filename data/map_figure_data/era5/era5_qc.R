library(tidyverse)
library(cowplot)

era5 <- read_csv("ERA_export_joe_QHI.csv") %>%
  mutate(month_long = case_when(month == 1 ~ "Jan",
                                month == 2 ~ "Feb",
                                month == 3 ~ "Mar",
                                month == 4 ~ "Apr",
                                month == 5 ~ "May",
                                month == 6 ~ "Jun",
                                month == 7 ~ "Jul",
                                month == 8 ~ "Aug",
                                month == 9 ~ "Sep",
                                month == 10 ~ "Oct",
                                month == 11 ~ "Nov",
                                month == 12 ~ "Dec",
                                )) %>%
  mutate(month_long = ordered(month_long)) %>%
  mutate(month_long = fct_reorder(month_long, month))
qa_plot <- ggplot(era5, aes(x = year, y= total_precipitation * 1000,
                   colour = as.character(month_long))) + geom_line() +
  labs(x = "", y = "ERA5 total monthly precipitation [mm]") +
  facet_wrap(vars(month_long), scales = "free_x") +
  theme_cowplot() +
  theme(legend.position = "none",
        strip.background = element_rect("NA"),
        axis.text.x = element_text(angle = 45, vjust = 0.5))
save_plot("era5_qc.png", qa_plot, base_height = 6)
