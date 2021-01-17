import pandas as pd
from Bio import Entrez
import GEOparse
import time
import json
import os
import logging

maxTriesCounter = 0


def setLogger(theLogger, fileName, logMode='a'):
       
    logsPath = r'C:\Users\shoshan1987\Desktop' 
    theLogger.handlers = []
    if not os.path.exists(logsPath):
        os.makedirs(logsPath, 0o755)
    logPath = os.path.join ( logsPath , fileName)
    rfh = logging.FileHandler(logPath, mode='w')
    formatter = logging.Formatter('%(asctime)s  [%(levelname)s]  %(filename)s::%(funcName)s  %(message)s')
    rfh.setFormatter(formatter)
    theLogger.addHandler(rfh)
    theLogger.setLevel(logging.DEBUG)
    theLogger.propagate = False
    pass

logger = logging.getLogger('logger')
setLogger(logger, "loggerShoshan_2.log", 'a')
logger.info("started running")
GSEIdNumber = "59847"



def GSEParsing(GSEIdNumber):
    #counter for not adding comma before first iteration    
    counter = 0
    GSEId = "GSE" + GSEIdNumber
    try:
        handleGSE = Entrez.esearch(db="gds", term="biopython", id=GSEId)
        recordGSE = Entrez.read(handleGSE)
        logger.info("Main data was loaded")
        print("Main data was loaded")
    except (RuntimeError, Exception) as e:
         logger.error(e, " Data wasn't loaded. trying again")
         print("There was a problem with the data loading process. Lets start from the begining")
         main(maxTriesCounter)
    
    GSEJson = '{"'+ GSEId+ '":[' + "\n"    
    
    try:
        for eachId in recordGSE['IdList']:
            counter = counter +1
            idFromIdList = str(eachId)
            handleId = Entrez.esummary(db="gds", id = idFromIdList)
            recordId = Entrez.read(handleId)
            Accession = recordId[0]['Accession']
            GPL = recordId[0]['GPL']
            suppFile = recordId[0]['suppFile']
            FTPLink = recordId[0]['FTPLink']
            jsonForLoadingToFile = converteDataToCSVFile(idFromIdList, Accession, GPL, suppFile, FTPLink)
            if counter != 1:
                GSEJson = GSEJson + ',' + "\n" + str(jsonForLoadingToFile) 
            else:
                 GSEJson = GSEJson + "\n" + str(jsonForLoadingToFile) 
            logger.info("Id data was loaded")
            print("Id data was loaded")     
    except Exception as e:
        logger.error(e, " Data wasn't loaded. trying again")
        print("There was a problem with the data loading process. Lets start from the begining")
        main(maxTriesCounter)
   
    GSEJson = GSEJson + ']'+ "\n" + '}'
    
    try:
        GSEJson = json.loads(eval(json.dumps(GSEJson)))
        logger.info("Created Json")
        print("Json was created")
    except Exception as e:
         logger.error(e, " Data didn't convert to json. trying again")
         print("There was a problem with the data converting to json process. Lets start from the begining")
         main(maxTriesCounter)   
    try:
        GSEJson = pd.DataFrame.from_dict(GSEJson, orient='columns')
        logger.info("Data Now is in pd df")
        
    except TypeError as e:
        logger.error(e, " Data didn't convert to json. trying again")
        print("There was a problem with the data converting to json process. Lets start from the begining")
        main(maxTriesCounter)   
    try:
        GSEJson.to_csv(GSEId +".csv", header = True)
        logger.info("Data uploaded to csv file")
        print("Data uploaded to csv file successfully")
    except AttributeError as e:
        logger.error(e, "Data didn't upload to csv file. trying again")
        print("There was a problem with the data cuploading to csv file process. Lets start from the begining")
        main(maxTriesCounter)   
    return(GSEId, idFromIdList, Accession, GPL, suppFile, FTPLink) 

def RnaSequencing(GSEIdw):
    #counter for not adding comma before first iteration
    counter =0
    deepSequencingDeatails = GSEIdw
    try:
        handleMain = Entrez.esearch(db="sra", term="seq", id=GSEIdw)
        recordMain = Entrez.read(handleMain)
        logger.info("Main data was loaded")
        print("Main data was loaded")
        
    except (RuntimeError, Exception) as e:
        logger.error(e, " Data wasn't loaded. trying again")
        print("There was a problem with the data loading process. Lets start from the begining")
        main(maxTriesCounter)
    try:
        for eachId in recordMain['IdList']:
            counter = counter +1
            idFromIdList = str(eachId)
            handleId = Entrez.esummary(db="sra", id = idFromIdList)
            recordId = Entrez.read(handleId)
            logger.info("Id data was loaded")
            print(("Id data was loaded"))
            if counter !=1:
                deepSequencingDeatails = deepSequencingDeatails + ',' + "\n" + str(recordId)
            else:
                deepSequencingDeatails = deepSequencingDeatails + "\n" + str(recordId)
    except Exception as e:
        logger.error(e, " Data wasn't loaded. trying again")
        print("There was a problem with the data loading process. Lets start from the begining")
        main(maxTriesCounter)
        
    try:
        FileName = GSEIdw + " DeepSequencingDeatails.csv"
        with open(FileName, "w") as file:
            file.write(deepSequencingDeatails)
        logger.info("Data uploaded to csv file")
        print("Data uploaded to csv file successfully")
        
    except AttributeError as e:
        logger.error(e, " Data didn't upload to csv file. trying again")
        print("There was a problem with the data cuploading to csv file process. Lets start from the begining")
        main(maxTriesCounter)       
    return


def converteDataToCSVFile(idFromIdList, Accession, GPL, suppFile, FTPLink):
    jsonForLoadingToFile=  "{" + '"id":"' + idFromIdList +'",'  + '"Accession":"' + Accession +'",' + '"GPL":"' + GPL + '",' + '"suppFile":"' + suppFile + '",' + '"FTPLink":"'  + FTPLink + '"}'
    return(jsonForLoadingToFile)


def main(maxTriesCounter):
    try:
        experimentSummaryOrDeepSequencing = input("Please type the desired output: Type 1 for Experiment Summary or Type 2 for Deep Sequencing ")
            
        if maxTriesCounter < 5:
            if (experimentSummaryOrDeepSequencing == "1") or (experimentSummaryOrDeepSequencing=="2"):
                logger.info("User chose" + experimentSummaryOrDeepSequencing)
                GSEOrPRJNA = input("Please type the desired GSE or PRJNA ")
                if "GSE" in GSEOrPRJNA or "PRJNA" in GSEOrPRJNA:
                    print("we are parsing your output, it may take a while")
                    if experimentSummaryOrDeepSequencing =="1":
                        try:
                            GSEParsing(GSEOrPRJNA)
                        except Exception as e:
                            logger.error(e, " Data wasn't loaded")
                        
                    else:
                        RnaSequencing(GSEOrPRJNA)
                else:
                    logger.info("Number of tries" + maxTriesCounter)
                    maxTriesCounter = maxTriesCounter + 1
                    print("The input must be in GSE or PRJNA type, Please restart from the begining")
                    time.sleep(2)
                    logger.info("Number of tries" + maxTriesCounter)
                    main(maxTriesCounter)        
            else:
                logger.info("Number of tries" + maxTriesCounter)
                maxTriesCounter = maxTriesCounter + 1
                print("The input must be 1 or 2")
                time.sleep(2)
                main(maxTriesCounter)
        else:
             print("you have reached to max tries, please call support")
             logger.info("Max tries")
         
    except Exception as e:
        print(e)
        logger.error(e, "Something went wrong")
    

if __name__ == "__main__":
    main(maxTriesCounter)

