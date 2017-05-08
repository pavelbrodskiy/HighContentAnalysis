@ECHO OFF
@set files=1
@set "activeDirectory=D:\Zartman Lab\Pavel\5.2.17 Cl8 Blebbistatin\"
@set "dirReset=%activeDirectory%reset.tif"
@set "dirHighFreq=%activeDirectory%PositionName_0001.tif"
@set "dirNewFolder=%activeDirectory%HighFrequency\"
:while1
ECHO.Checking for output file
IF EXIST "%dirReset%" (
    ECHO.File Exists
    copy "%dirHighFreq%" "%dirNewFolder%HighFreq_%files%.tif"
    @set /a files=files+1
    del /f "%activeDirectory%template2.STG"
    copy "%activeDirectory%template.STG" "%activeDirectory%template2.STG"
    more +5 "%activeDirectory%positions2.STG" >> "%activeDirectory%template2.STG"
    del /f "%activeDirectory%positions2.STG"
    move "%activeDirectory%template2.STG" "%activeDirectory%positions2.STG"
    del /f "%dirReset%"
) else (
    ECHO.File Doesn't Exist
)
ECHO.Iteration number %files%
TIMEOUT 5
goto :while1
ECHO ON