
# Amanda

I've added the individual linear regression lines in the background. However, with the two colors/linetypes it's very difficult to distinguish the main two lines for the groupings. This is an image that would go into Figure 1. Also it's for MI/15 instead of CA/09 so there about 200 less individuals (400 less lines considering pre/post titer) present on these plots compared to what would be later for CA/09. Also, with the change in the HAI virus panel the domains of the linear regressions are different. Although it provides more transparency in that aspect it does add to the confusion. (Fig 1)

For the bounding lines, I've calculated the linear regression for each person, pre and post vaccination. Below are the averages and sd of the slopes and the intercepts. My question is how do I incorporate the intercept? If it was only the slope I understand that it would just be above and below, but the intercept is not fixed. Would there be 4 potential bounding lines, or am I thinking about this incorrectly? (Fig 1 alt)

I'll include the bounding lines for the titer = slope_upperSD * distance + intercept_upperSD (purple line) and opposite for lower bound (green line) since they include the other combinations (slope_lwr + intercept_upr (pink line)  (Fig 2)

Lo siento about the spam of emails, but we were wrong in the previous meeting. I didn't compare the overall linear regression results to averaging the independent linear regressions together. I looked it up in my notes and it turns out what we were remembering was the average AUC of the individuals being equal to the AUC of the linear regression. 

I'm mentioning this because a problem arises when we take the average of the linear regressions. The linear regression and the averaged Linear regression do not match. If the ranges of the data were the same it wouldn't be a problem but the data that does not have as large of a range have a steeper slope and are pulling the averaged linear regression slope results to be steeper. This isn't matched in the linear regression with all of the data since it doesn't extrapolate for those points. (Fig 3)

# Zane
pdf doc

# Andreas

My thoughts/comments:

* Is the reason that we have different amounts of data spanning the x-axis differently because of different seasons? For a given season/year, every person should have been tested on the same panel of strains, correct? So if we did a per-season analysis, would this unbalanced data go away? (I’m not saying we should switch to that analysis only, just wondering how this unbalance comes about). And also, like Zane, I’m wondering about this systematic steeper slope for those with less data. I can’t fully understand just now why that happens. And some of the lines go systematically up as we go further. That seems curious too. Might be worth poking at the data for a few of those individuals to see what’s happening there.

* Do we still have that AUC comparison? Agree with Zane that unclear why they should give the same results. Just curious to take another look to see how that happens.

* Conceptually, I like the hierarchical/partial pooling approach the best to get some overall estimate and uncertainty around it. I was just trying to stay away from it since Amanda I don’t think you’ve done any hierarchical/mixed-effects modeling and I didn’t want to throw some new methodology at you. But maybe if the 2 of you work on that, Zane could write down/implement the model? (Of course I prefer if it were Bayesian, but I don’t want to change things around too much, so if we want to do ML for now, ok with me.)

* Is everything now in the Onedrive folder and can everyone access the latest? I might want to play around a bit with it myself, and I want to make sure we are now all working in the same space.

* Of course, happy to have another call/conversation if useful. Or we keep using email. Either is fine with me.

# Amanda

* Yes, since the graphs include different seasons the panels were different which affects the maximum distance that some individuals have. The unbalanced would go away by looking at it only for season/year. 
* The AUC comparison is in the exploratory folder I believe. I'm pretty sure that is what I was remembering. (~ALS033 - Antigenic Distance/Dropbox Contents/code/3_exploratory/7_Unaggregated AUC_of_individuals.Rmd )
* I would be ok with implementing a different technique if it is more appropriate than the method we're currently doing. I would need to defer the work to Zane however.
* Everything is in the Onedrive, and should be syncing appropriately.

<!-- END OF FILE -->
