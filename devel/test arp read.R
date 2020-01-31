rm(list = ls())
library(strataG)
library(tidyverse)

x1 <- arlequinRead("good_ne_test_1_1.arp")
x2 <- arlequinRead("problem_ne_test_1_2.arp")

g1 <- arp2gtypes(x1)
g2 <- arp2gtypes(x2)
