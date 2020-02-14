# findsaccades
The manual saccade marking is done through findsaccades and extract saccades code. Once that is complete you can move onto outputting it to an excel sheet. 

Excel sheet generation with rethresholding microsaccade start and end periods:

Run gothroughalldatafixfast first using the velocity figures and the raw outputs from TSLO motion extraction
Run getalldataforsheet selecting the output files from the previous code to generate an excel sheet with data for each individual microsaccade.
Run getaveragesexcelsheet to generate a sheet where each row is the average of all traces for each patient. 
