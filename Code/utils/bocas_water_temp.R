
library(tidyverse)
library(lubridate)

read_csv('../../Data/bocas_temp/bocas_tower_wt_elect.csv') %>%
  select(datetime, wt, chk_note) %>%
  mutate(datetime = mdy_hms(datetime),
         year = year(datetime))
