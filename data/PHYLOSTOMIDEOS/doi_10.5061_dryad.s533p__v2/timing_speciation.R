# R script to quantify the number and timing of speciation events between extant
# sister species for the Quaternary and pre-Quaternary (Neogene + Paleogene)
# periods on a posterior sample of time-calibrated phylogenies. The code can be
# modified to accommodate older periods if necessary.

library(ape)

# Load the file with the posterior distribution of time-calibrated trees. In
# this example we named the file as "my_trees.nex".

trees <- read.nexus("my_trees.nex")

# Create three lists.
# bal stores the number of descendants (i.e. tips) on each of the daughter-
# branches for each node.
# brt stores the distance from each node to the tips.
# tip.br stores the branch lengths of the tips.

bal <- list()
brt <- list()
tip.br <- list()

# Create three lists and a matrix to store the number of speciation events of
# extant species for the Quaternary (Qua), the Neogene (Neo) and the Paleogene
# (Pal).

Qua <- list()
Neo <- list()
Pal <- list()

counts <- data.frame(matrix(, length(trees), 3, dimnames = list(c(1:length(trees)), c("Qua", "Neo", "Pal"))))

# Calculate the number of speciation events.

for (i in 1:length(trees)){
  bal[[i]] <- balance(trees[[i]])
  brt[[i]] <- branching.times(trees[[i]])
  tip.br[[i]] <- list()
  for (j in 1:trees[[i]]$Nnode){
    if (bal[[i]][j, 1] == bal[[i]][j, 2]){
      if (bal[[i]][j, 1] == 1){
        tip.br[[i]][[j]] = brt[[i]][[j]]
      }
    }
  }
  tip.br[[i]] <- sapply(tip.br[[i]], function(x){as.numeric(x)[1]})
  tip.br[[i]] <- tip.br[[i]][!is.na(tip.br[[i]])]
  Qua[[i]] <- sum(tip.br[[i]] < 2.588)
  Neo[[i]] <- sum(tip.br[[i]] >= 2.588 & tip.br[[i]] < 23.03)
  Pal[[i]] <- sum(tip.br[[i]] >= 23.03 & tip.br[[i]] < 66)
  counts[i, 1] <- Qua[[i]]
  counts[i, 2] <- Neo[[i]]
  counts[i, 3] <- Pal[[i]]
}

# Since we are interested in Quaternary and pre-Quaternary events, we sum the
# number of events of the Neogene and Paleogene (i.e. pre-Quaternary).

counts$pre.Qua = counts$Neo+counts$Pal

write.table(counts, "speciation_events.txt", sep = "\t", row.names = F)

# Calculate the timing of speciation events and store the results in a table.

timing <- as.numeric(unlist(tip.br))
timing.Qua <- timing[which(timing < 2.588)]
timing.Neo <- timing[which(timing >= 2.588 & timing < 23.03)]
timing.Pal <- timing[which(timing >= 23.03)]
timing.pre.Qua <- c(timing.Neo, timing.Pal)

timing.data <- data.frame(matrix(, length(timing), 2, dimnames = list(c(1:length(timing)), c("Period", "Timing"))))
timing.data[, 1] <- rep(c("Quaternary", "pre Quaternary"), c(length(timing.Qua), length(timing.pre.Qua)))
timing.data[1:length(timing.Qua), 2] <- timing.Qua
timing.data[(length(timing.Qua)+1):length(timing), 2] <- timing.pre.Qua

write.table(timing.data, "timing_of_events.txt", sep = "\t", row.names = F)

# Plot the results for the number of speciation events.

library(ggplot2)

new.table <- data.frame(matrix(, 2*length(trees), 2, dimnames = list(c(1:(2*length(trees))), c("Period", "Events"))))
new.table[, 1] <- rep(c("Quaternary", "pre Quaternary"), c(length(trees), length(trees)))
new.table[1:length(trees), 2] <- counts$Qua
new.table[(length(trees)+1):(2*length(trees)), 2] <- counts$pre.Qua

p1 <- ggplot(new.table, aes(x = Events, fill = Period)) + geom_histogram(position = "identity", alpha = 0.75, binwidth = 1)
p1 <- p1 + theme_bw()
p1 <- p1 + scale_fill_manual(values = c("#F5A275", "#99C7EC"))
p1 <- p1 + xlab("Number of speciation events (extant species)") + ylab("Counts")
print(p1)
ggsave("number_of_events.pdf", width = 16, height = 14, units = "cm")

# The End
