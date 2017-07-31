#!/usr/bin/env python

# A python program to identify the amount of games that you can influence assuming that both teams are
# equally matched skillwise. Since you have x (11 for Overwatch) other players playing in the same game,
# whether you win or not is not completely determined by how well you play, it is heavily dependent on
# how well your teamates play. This program determines the amount of games where you can change the
# outcome of a game by playing well or badly. This does not take into account draws.

import matplotlib.pyplot as plt
import random
import numpy as np
import scipy.stats
from scipy.stats import truncnorm

# These are all the user defined variables, the number of Monte Carlo statistics, the distributions we
# wish to consider, and the number of other players in the game e.g. 11 for Overwatch, 9 for Dota
M=100000
distribution_names=["discrete","uniform","truncated normal distribution with scale 0.25", \
        "truncated normal distribution with scale 0.5","truncated normal distribution with scale 1", \
        "truncated normal distibution with scale 2", "truncated normal distribution with scale 4"]
distributions=["random.choice([-1,0,1])","random.uniform(-1,1)","truncnorm.rvs(-4,4,scale=0.25)", \
        "truncnorm.rvs(-2,2,scale=0.5)","truncnorm.rvs(-1,1,scale=1)", \
        "truncnorm.rvs(-0.5,0.5,scale=2)","truncnorm.rvs(-0.25,0.25,scale=4)"] 
number_of_other_players=11 # 11 for Overwatch, 9 for Dota etc

scores=[0]*len(distributions)
definite_losses=[0]*len(distributions)
definite_wins=[0]*len(distributions)
games_you_can_influence=[0]*len(distributions)
monte_carlo_distribution=[[0]*M for i in range(len(distributions))] # has dimensions MCD[len(dist)][M]

print "Using",M,"as the number of Monte Carlo statistics."
print "Using",number_of_other_players+1,"as the total number of players."
print "Using",len(distributions),"different distributions."
print "Printing progress of the Monte Carlo simulation:\n"

for n in range(M):
    for other_players in range(number_of_other_players):
        for count, dist in enumerate(distributions):
            monte_carlo_statistic=eval(dist)
            monte_carlo_distribution[count][n]=monte_carlo_statistic
            scores[count]=scores[count]+monte_carlo_statistic
    for count, score in enumerate(scores):
        if 1 < score <= number_of_other_players:
            definite_wins[count]=definite_wins[count]+1
        elif (-1)*number_of_other_players <= score < -1:
            definite_losses[count]=definite_losses[count]+1
        elif -1 <= score <= 1:
            games_you_can_influence[count]=games_you_can_influence[count]+1
        else:
            exit("Score too big or too small.")
    scores=[0]*len(distributions)
    # Add a progress indicator
    if n!=0:
        if n%(M/100)==0:
            print (100*n)/M,"\b% done."

print ""
for count, dist in enumerate(distribution_names):
    print dist.upper()
    print "Percentage of games that are definite wins is",100*float(definite_wins[count])/M,"\b%"
    print "Percentage of games that are definite losses is",100*float(definite_losses[count])/M,"\b%"
    print "Percentage of games that you can influence is",100*float(games_you_can_influence[count])/M,"\b%\n"
    print "Now assuming perfect play (you always score a 1):"
    print "Percentage of games that you will win is",100*(float(definite_wins[count])+ \
        games_you_can_influence[count])/M,"\b%"
    print "Percentage of games that you will lose is",100*float(definite_losses[count])/M,"\b%\n"
    # test values for the bw_method option ('None' is the default value)
    bw_values =  [0.1]
    # generate a list of kde estimators for each bw
    kde = [scipy.stats.gaussian_kde(monte_carlo_distribution[count],bw_method=bw) for bw in bw_values]
    # plot (normalized) histogram of the data
    plt.hist(monte_carlo_distribution[count], 50, normed=1, facecolor='green', alpha=0.5);
    # plot density estimates
    t_range = np.linspace(-2,8,200)
    for i, bw in enumerate(bw_values):
        plt.plot(t_range,kde[i](t_range),lw=2)
    plt.xlim(-1,6)
    plt.legend(loc='best')
    axes = plt.gca()
    axes.set_xlim([-1.2,1.2])
    axes.set_ylim([0,2])
    filename='figure_' + str(number_of_other_players) + '_' + str(M)  + '_' + str(count) + '.png'
    plt.savefig(filename)
    print "Saved histogram of Monte Carlo data to",filename,"\n\n"
    plt.clf()
