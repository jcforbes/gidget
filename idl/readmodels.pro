
FUNCTION readmodels,nameList
        model=readOutput(nameList[0])
        modelList=ptrarr(n_elements(nameList), /allocate_heap)
        *(modelList[0])=model
        FOR i=1, n_elements(nameList)-1 DO *(modelList[i])=readOutput(nameList[i])
        RETURN,modelList
END
