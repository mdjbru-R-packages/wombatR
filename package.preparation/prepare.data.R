# data files are taken from
# http://www.wildanimalmodels.org/tiki-index.php?page=The+ecologists+guide+to+the+animal+model#WOMBAT
# which is the wiki corresponding to the paper
# Wilson, Alastair J., Denis Reale, Michelle N. Clements, Michael M. Morrissey,
# Erik Postma, Craig A. Walling, Loeske E. B. Kruuk, and Daniel H. Nussey. “An
# Ecologist’s Guide to the Animal Model.” Journal of Animal Ecology 79, no. 1
# (January 2010): 13–26. doi:10.1111/j.1365-2656.2009.01639.x.



ped = read.table("gryphon.ped_WOMBAT")

names(ped) = c("animal", "sire", "dam")

data = read.table("gryphon_uni.dat")

names(data) = c("animal", "dam", "birth.year", "sex", "birth.weight")

data = data[, c(1, 3:5)]

gryphon = merge(x = ped, y = data, by = "animal", all = T)

save(gryphon, file = "gryphon.rda")

