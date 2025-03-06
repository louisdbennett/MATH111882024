rm(list = ls(all = TRUE))
# Load data
unis <- read.csv("finalbis.csv")


# Ethnic Diversity Index Computation: https://www.ed-data.org/article/Ethnic-Diversity-Index

# Add Ethnic Diversity Index
c1 = 100 
c2 = -100*sqrt(5*(5-1))/(5-1)
unis$ethnic_diversity_index = c1 + c2*sqrt((unis$White.ethnic.group - 0.2)^2 + (unis$Black.ethnic.group - 0.2)^2
                                            + (unis$Asian.ethnic.group - 0.2)^2 + (unis$Mixed.ethnic.group - 0.2)^2
                                            + (unis$Other.ethnic.group - 0.2)^2)

# Add Russell Group Uni factor
unis$russell = FALSE
russell_codes = c("B32", "B78", "C05", "C15", "D86", "E56", "E84", "G28", "I50", 
                  "K60", "L23", "L41", "L72", "M20", "N21", "N84", "O33", "Q50", "Q75",
                  "S18", "S27", "U80", "W20", "Y50")
unis[unis$INSTITUTION_CODE %in% russell_codes, ]$russell = TRUE

unis

# Task: Predicting graduate employability and the most predictive features
# Part 1: Basic Linear Regression Model
lm1 <- lm(career_after_15_month ~ Women + russell + ethnic_diversity_index + POLAR4.Q5, data = unis)
summary(lm1)