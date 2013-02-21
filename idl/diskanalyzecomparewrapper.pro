PRO diskAnalyzeCompareWrapper,nameList,sv,cmt,leg
	IF(n_elements(nameList) EQ 0 || n_elements(sv) EQ 0 || n_elements(cmt) EQ 0 || n_elements(leg) EQ 0) THEN BEGIN
		message,"Usage: diskAnalyzeCompare,['<name1>','<name2>',...],sv,comment,['<leg1>','<leg2>',...] -- sv=0->don't save sv=1->save as mpg sv=2->save png images to be constructed into a movie by ffmpeg, sv=3->same as 2 but don't plot in an X window."
	ENDIF

	modelList = ReadModels(nameList)
	
	diskAnalyzeCompare,modelList,sv,cmt,leg
END

