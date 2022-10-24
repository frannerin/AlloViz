# Major problems
- There are no tests (normal or regression tests)
- There is no documentation nor list about analysis methods other than btw. 
- Analysis method raw not documented
- The API makes no sense to me. Why storing multiple calculations in the same object? 
- Having many ways to do the same computation has no logic and is error-prone. THere are many ways to filter and I can't even find them any more.
- Filtering *methods* should start lower case
- Detect GPCRs when they don't come from GPCRmd
- Data caching is confusing. What if one changes the underlying trajectory and calls compute() again? Is the "data" dir reused?

# Packaging problems
- Make a list of imports without versions, or it can't be used
- Removing files from packages confuses the history, impedes updates from upstream, and saves no space. May also be seen as unnecessary tweaks to the original work.
- By all means keep the packages as they are, if possible. Tweak the path instead of making relative imports
- Relies on vmd-python which violates VMD license


# Minor
- Lazy import via exec is weird and delays discovery of installation problems

# Resolved
- getcontacts relies on vmd-python which seems outdated. I made a patch to have it work on new macs
