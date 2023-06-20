library(devtools);library(tidyverse);library(lubridate)

use_build_ignore("/ignore_things/ignored.R")


#### Create dummy dataset ####

read_csv("/Users/jonas/OneDrive - Syddansk Universitet/Lyngsø/Forze/CH4_corrected.csv") %>%
  filter(wp == 49) -> data

data %>%
  filter(PumpCycle %in% c(9,13)) %>%
  group_by(PumpCycle) %>%
  mutate(PumpCycle = if_else(PumpCycle == 9, 1,2),
         time = datetime-min(datetime),
         datetime = ymd_hms("2021-09-28 03:11:35") + time + (PumpCycle-1)*3600,
         pred_CH4 = pred_CH4-15,
         station = 1,
         sensor = 1) %>%
  select(datetime,`RH%`,tempC,CH4smV,K33_RH,K33_Temp,K33_CO2,SampleNumber,PumpCycle,pred_CH4,station,sensor) -> DIY_sensor_data

DIY_sensor_data %>%
  ggplot(aes(datetime,pred_CH4, group = PumpCycle)) +
  geom_line()

use_data(DIY_sensor_data, overwrite = T)
