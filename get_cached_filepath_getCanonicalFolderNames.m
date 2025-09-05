function folderName = get_cached_filepath_getCanonicalFolderNames(folderName)
folderName = regexprep(folderName, '^\\\\ucibciserver\.d\.zind\.org\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
folderName = regexprep(folderName, '^\\\\ucibciserver\.myqnapcloud\.com\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
folderName = regexprep(folderName, '^\\\\128\.200\.248\.115\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
folderName = regexprep(folderName, '^\\\\10\.141\.150\.188\\', '\\\\ucibciserver\.stemcell\.uci\.edu\\');
