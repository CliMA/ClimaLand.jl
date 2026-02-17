update_photosythesis!



is_c3
weighted average of c3 and c4

c3 + c4 = 1

c3_mixed keyword argument


Tests might be testing linear interpolation? I am not sure...

Steps / approach for doing this PR
- Investigate whether the values are 0 and 1 or a continuous spectrum
- is_c3 should be replaced with fractional_c3
- Patterns of the form is_c3 > 0.5 ? blah1 : blah2 should be replaced with a weighted average
- result depending on is_c3 should also be replaced with weighted average too (not all functions are linear though!)
	- This part was a little confusing to me, so might need to ask about this
	- The example was that if J depends on c3 and c4 and J is used elsewhere, you can't just pass in J to compute idk. You need to
	compute idk by doing a weighted average of J_c3 and J_c4 (both computed separately with c3 and c4).
- Backward compatibility is maintained because we are doing weighted averages
- Should take about four days?
