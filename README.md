# Reproducible research: version control and R

## Questions 1,2 and 3:
https://github.com/1076644/logistic_growth_computer_week3/tree/main

## Q4) (30 points) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:
### a. A script for simulating a random_walk is provided in the question-4-code folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points)

If we execute to code given in the question-4-code folder, we produce two plots side by side simulating a random walk for two independent iterations. 500 steps are performed in each plot, with the gradient of the colour of blue dictating how many of the 500 steps have been taken. This allows us to see how the random motion changes over time, with earlier steps being shown in a darker blue which progressively becomes lighter. the sample unit for each step is set to 0.25, with a random angle taken from 0 to 2π radians. Combining these metrics over 500 iterations generates x and y coordinates in both the plots. reproducing this code over and over again however produces random plots that are not repeatable, as a different set of random angles are used whihc produce different trajectories. This leads to issues if wanting to use thos code over and over again, or by giving it to someone else to see if they can copme to the same results. We do see both plots beginning at the origin, but there is no observable relationship between the graphs, which is a common characteritics of Brownian motion. 

![image](https://github.com/user-attachments/assets/56a0f454-59eb-4b28-a626-19311c334759)

- As we can see, we produce two independent iterations of the code, and if re-run, these graphs would not be reproduced, which gives us many issues when wanting to share our code or come back and change things later.

### b. Investigate the term random seeds. What is a random seed and how does it work? (5 points)

A random seed is a number or vector that is used to intialise random number generators. The random number generator is completely determined by the seed with the random seed allowing other researchers to replicate the results exactly, which is highly desirable in the scientific world. This is especially crucial when replicating random processes such as Brownian motion, whihc is being used as an example here.

### c. Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked reproducible-research_homework repo. (10 points)

The original script for the Brownian motion can be shown below:

``` r
#install.packages("ggplot2")
#install.packages("gridExtra")

library(ggplot2)
library(gridExtra)

random_walk  <- function (n_steps) {
  
  df <- data.frame(x = rep(NA, n_steps), y = rep(NA, n_steps), time = 1:n_steps)
  
  df[1,] <- c(0,0,1)
  
  for (i in 2:n_steps) {
    
    h <- 0.25
    
    angle <- runif(1, min = 0, max = 2*pi)
    
    df[i,1] <- df[i-1,1] + cos(angle)*h
    
    df[i,2] <- df[i-1,2] + sin(angle)*h
    
    df[i,3] <- i
    
  }
  
  return(df)
  
}

data1 <- random_walk(500)

plot1 <- ggplot(aes(x = x, y = y), data = data1) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

data2 <- random_walk(500)

plot2 <- ggplot(aes(x = x, y = y), data = data2) +
  
  geom_path(aes(colour = time)) +
  
  theme_bw() +
  
  xlab("x-coordinate") +
  
  ylab("y-coordinate")

grid.arrange(plot1, plot2, ncol=2)
```
- This code prodcues the random graphs explained in Q4)a.

The edited script can be observed in the file "Edited Script.Rmd" in the `question-4-code` file in the `reproducible-research_homework` repo. Evidence of this script now producing reproducible results is shown in the figure below where seperate data sets used for Brownian motion produce the same graph due to the addition of a seed number. 

#### Reprducible Figure 
![image](https://github.com/user-attachments/assets/1f16032b-1e2f-4c43-a0f4-df5c3bf78699)

### d. Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the README.md of the fork). (5 points)

Below shows a screenshot of my commit history to the Edited code script and the changes made from the original `question-4-code` provided to ensure that this code is now reproducible.

#### Commit history illustrating the changes made to the code
![image](https://github.com/user-attachments/assets/fb3975fb-cbf1-4c26-ba50-dc52c89faf7f)

## Q5) (30 points) In 2014, Cui, Schlub and Holmes published an article in the Journal of Virology (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \alpha L^{\beta}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

### a. Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)

The code below demonstrates how we can obtain the data, and then see how many rows and columns the table has 
 ```r
DNAdata <- read.csv("Cui_etal2014.csv")

nrow(DNAdata) # This gives us the number of rows in the data set
ncol(DNAdata) # Thsi gives us the number of the columns in the data set 
```
- The output from this gives us **33 rows** and **13 columns**

### b. What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points)

For this data set, a natural log transformation can be done, which can change the allometric equation into one which resembles the form **$`y = ax + b`$**. 

```math
V = \alpha L^{\beta}
```
This can be transformed into:
```math
ln(V) = ln(\alpha) + {\beta}ln(L)
```

This then resembles a linear model which can then be fitted to this data. This can be done by using the following code below:

```r
DNAdata$Virion_volume_log <- log(DNAdata$Virion.volume..nm.nm.nm.)  # Log-transformed Virion volume
DNAdata$Genome_length_log <- log(DNAdata$Genome.length..kb.)  # Log-transformed Genome length

head(DNAdata) #To check this has formatted properly
```
Applying the linear model to the data then looks like this:

```r
linear_regression <- lm(Virion_volume_log ~ Genome_length_log, DNAdata)
summary(model1)
```
The output from this looks as follows:

![image](https://github.com/user-attachments/assets/730b3a84-ab2d-4996-95d8-d17202ccedd5)

### c. Find the exponent ($\beta$) and scaling factor ($\alpha$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points)

From the table above of our linear model, we can find the exponent ($\beta$) and the scaling factor ($\alpha$), as well as the p values for these estimates to see how good our transformation is for obtaining these key values. 

- Our scaling factor according to the linear model is 7.0748, which is hswon through the intercept. However, this is the ln value, so to get the true value, we need to do **$`e^{7.0748}`$** which gives us a value of **1181.807116**
## Instructions

The homework for this Computer skills practical is divided into 5 questions for a total of 100 points. First, fork this repo and make sure your fork is made **Public** for marking. Answers should be added to the # INSERT ANSWERS HERE # section above in the **README.md** file of your forked repository.

Questions 1, 2 and 3 should be answered in the **README.md** file of the `logistic_growth` repo that you forked during the practical. To answer those questions here, simply include a link to your logistic_growth repo.

**Submission**: Please submit a single **PDF** file with your candidate number (and no other identifying information), and a link to your fork of the `reproducible-research_homework` repo with the completed answers (also make sure that your username has been anonymised). All answers should be on the `main` branch.

## Assignment questions 

1) (**10 points**) Annotate the **README.md** file in your `logistic_growth` repo with more detailed information about the analysis. Add a section on the results and include the estimates for $N_0$, $r$ and $K$ (mention which *.csv file you used).
   
2) (**10 points**) Use your estimates of $N_0$ and $r$ to calculate the population size at $t$ = 4980 min, assuming that the population grows exponentially. How does it compare to the population size predicted under logistic growth? 

3) (**20 points**) Add an R script to your repository that makes a graph comparing the exponential and logistic growth curves (using the same parameter estimates you found). Upload this graph to your repo and include it in the **README.md** file so it can be viewed in the repo homepage.
   
4) (**30 points**) Sometimes we are interested in modelling a process that involves randomness. A good example is Brownian motion. We will explore how to simulate a random process in a way that it is reproducible:

   a) A script for simulating a random_walk is provided in the `question-4-code` folder of this repo. Execute the code to produce the paths of two random walks. What do you observe? (10 points) \
   b) Investigate the term **random seeds**. What is a random seed and how does it work? (5 points) \
   c) Edit the script to make a reproducible simulation of Brownian motion. Commit the file and push it to your forked `reproducible-research_homework` repo. (10 points) \
   d) Go to your commit history and click on the latest commit. Show the edit you made to the code in the comparison view (add this image to the **README.md** of the fork). (5 points) 

5) (**30 points**) In 2014, Cui, Schlub and Holmes published an article in the *Journal of Virology* (doi: https://doi.org/10.1128/jvi.00362-14) showing that the size of viral particles, more specifically their volume, could be predicted from their genome size (length). They found that this relationship can be modelled using an allometric equation of the form **$`V = \alpha L^{\beta}`$**, where $`V`$ is the virion volume in nm<sup>3</sup> and $`L`$ is the genome length in nucleotides.

   a) Import the data for double-stranded DNA (dsDNA) viruses taken from the Supplementary Materials of the original paper into Posit Cloud (the csv file is in the `question-5-data` folder). How many rows and columns does the table have? (3 points)\
   b) What transformation can you use to fit a linear model to the data? Apply the transformation. (3 points) \
   c) Find the exponent ($\beta$) and scaling factor ($\alpha$) of the allometric law for dsDNA viruses and write the p-values from the model you obtained, are they statistically significant? Compare the values you found to those shown in **Table 2** of the paper, did you find the same values? (10 points) \
   d) Write the code to reproduce the figure shown below. (10 points) 

  <p align="center">
     <img src="https://github.com/josegabrielnb/reproducible-research_homework/blob/main/question-5-data/allometric_scaling.png" width="600" height="500">
  </p>

  e) What is the estimated volume of a 300 kb dsDNA virus? (4 points) 
