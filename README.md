# Presidential_election_forecast_2016
A Bayesian analysis of polling data and historical voting records to try to predict the outcome of the 2016 U.S. Presidential race.  
Carl Ehrett, November 7 2016

## Introduction

In this project I present a Bayesian model for predicting the outcome of the 2016 presidential election based on a combination of historical election results (from the previous 10 elections) and polling data (from Google Consumer Surveys). 

## Model Details 
The model views the election in each state $i$ as an outcome of a bernoulli $p_i$ trial. A set of (roughly) weekly polls in each state is also treated each as an outcome of a $p_i$ trial. To accommodate shifts in $p_i$ throughout the course of the campaign, these polls are included by way of a power prior with exponent $a_0$. The use of the power prior allows polls closer to the election to weigh more heavily than earlier polls. Electoral history in each state is used to form an informative prior on $p_i$. Using Metropolis-Hastings-within-Gibbs Markov-chain Monte Carlo sampling, values for $p_i$ and $a_0$ are then estimated. The values for $p_i$ are used to draw an outcome for the election.

## Conclusion
This project was completed in early November 2016, prior to the election. The forecast gives Hillary Rodham Clinton a 54% chance of becoming President, and Donald Trump a 46% chance.

## Full details
For full details of the model, the data sources, the methodology used, and the results of the analysis, please consult the report [here](./pres_elec_forecast.pdf "Full report"). The report includes a breakdown by state, showing for each state the estimated probability of Clinton victory. The report also indicates which states should be considered swing states.

## Repository contents
The contents of this github repository are as follows:  
[main.r](./main.r "Main R script"): Main R script, which performs the Bayesian analysis.  
[all_polls.csv](./all_poll.csv "Google Consumer Surveys data"): Collection of 13 Google Consumer Survey polls from Aug. 10 to Nov. 1 2016.  
[modern_results_by_state.csv](./modern_results_by_state.csv "Historical election data"): Spreadsheet of presidential electoral results by state and election year, from 1976 to 2012.  

## To run the model
To run the software with respect to the 2016 election using the included data from the Google Consumer Survey polls, simply run the R script [main.r](./main.r) in a working directory containing the spreadsheets [all_polls.csv](./all_poll.csv "Google Consumer Surveys data") and [modern_results_by_state.csv](./modern_results_by_state.csv "Historical election data"). 

To run the software with respect to another Presidential election, simply replace the two spreadsheets with spreadsheets containing polling data and historical data relevant to the election in question. The historical data should take the form (like modern_results_by_state.csv) of a sheet in which each row corresponds to a state (in alphabetical order), each column corresponds to an election year (in chronological order), and each entry contains a 1 for Democratic victory in that state and year, and a 0 for Republican victory. (No third party has won electoral votes in recent history; if this changes, substantive revisions will be required for this model to be used.) The polling data should conform to the format of the Google Survey polls (for which, see either [the full report](./pres_elec_forecast.pdf "Full report") or else directly examine the spreadsheet [all_polls.csv](./all_poll.csv "Google Consumer Surveys data")). Alternatively, polling data in a different format could be used, with modifications made to [main.r](./main.r "Main R script") in order to accommodate the new format.
