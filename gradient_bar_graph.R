library(tidyverse)
library(dplyr)
df <- read.table("B_mol_f.txt", sep = "\t", header = T)
library(lubridate)
View(df)

  ggplot(df, aes(reorder(GO_molecular_function_complete, Counts), Counts)) + 
  geom_col(aes(fill = P_value)) + 
  scale_fill_gradient2(mid  = "red", 
                      high = "blue", 
                       midpoint = mean(df$P_value)) + 
  coord_flip() + 
  labs(x = "molecular_function")
