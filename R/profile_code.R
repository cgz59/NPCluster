function() {
	prof <- profileExample(n = 50, p = 250, n.burn = 100, n.reps = 200)
	rows <- grep("elementwise_DP.functions.R", row.names(prof))
	
	
	prof <- prof[order(expand_numbers(row.names(prof))),]
	
	
}


expand_numbers <- function(names) {
	names <- sub("R#(\\d)$", "R#000\\1", names)
	names <- sub("R#(\\d\\d)$", "R#00\\1", names)
	names <- sub("R#(\\d\\d\\d)$", "R#0\\1", names)
	names
}