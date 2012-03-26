;; a function to convert from lookback time to redshift using available data
;; i.e. a model object stored in the common block modelBlock
FUNCTION conv_axis,axis,index,value
        COMMON modelBlock
        ti=nearest("time",1,theModel,value)
        z=theModel.evArray[9,ti]
        RETURN,string(FORMAT='(F0.1)',abs(z))
END