top-makers are the 250 markers with the lowest correlation matrix row sum
(i.e. I calculated a correlation matrix of g, applied an abs() on it and then summed the rows/cols, as a simple and quick metric of how much information they add (of course not ideal or technically sound))
Then I added the ones selected in the GA.
For top_markers2 I added 1000 random ones on top of that (and removed duplicates)
