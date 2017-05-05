fchoose <- file.choose()
slash <- gregexpr("/", fchoose)[[1]]
setwd(substring(fchoose, 1, slash[length(slash)]))

setwd("/media/jeff/Médiathèque/Audio/Chiroptera/00_Backup_Recordings/0_Analyse_OK/20160526_LIFE_Rochefort")
#setwd("/home/jf/R/SOUNDS/batmetrics")
#rm(list=ls())
path <- getwd()
path.stspt <- strsplit(path,"/")
fol.name <- path.stspt[[1]][length(path.stspt[[1]])]
## Load wav file
list.files(pattern = "wav",ignore.case=T)

require(seewave)
require(tuneR)
#### Change filename and comment
filename <- list.files(pattern = "wav",ignore.case=T)[2]
PI <- "TrainMYODAU_20170505"
smpl <- readWave(filename) # tuneR::readWave()
#wave <- smpl
#tr <- timer.abs(smpl@right, f=44100, threshold=10, msmooth = c(8820,50)) #resolution:2sec
#(Nchunk <- length(tr$s.start))
#names.chk <- paste0("Chunk_",1:Nchunk)
### make equal length of start and end vectors by adding -if unequal- total length
#if(Nchunk > length(tr$s.end)){
#  tr$s.end[Nchunk] <- length(smpl@right) / 44100}
#chnk <- list()
#rangeTE <- c()
### Loop detection of signals
#for(j in 1:Nchunk){
#  chnk[[j]] <- cutw(smpl@right, f=44100,
#                    from= tr$s.start[j], to=tr$s.end[j])
#  mainTE.rg <- range(chnk[[j]])
#  rangeTE[j] <- mainTE.rg[2] - mainTE.rg[1]
#}
#data.frame(CHUNK=names.chk,
#           RANGE=rangeTE,
#           Start.sec=round(tr$s.start,2),
#           Start.MinSec=paste0(tr$s.start %/% 60,"'",round(tr$s.start %% 60,1),"''"))

# Detect noise level
tr.noise <- timer(smpl@left, f=44100, threshold=5,msmooth=c(220.5, 90), plot=T)
length(tr.noise$s.start)
Mx <- numeric()
for(i in 1:30){
  Mx[i] <- median(env(cutw(channel(smpl, "left"),from=tr.noise$s.start[i] - 0.06
                           , to=tr.noise$s.start[i] - 0.04, output = "Wave"), plot=F, envt="abs"))}
LevelsAll <- fivenum(env(channel(smpl, "left"), plot=F, envt="abs"))
Mx
LevelsAll
(noise <- fivenum(Mx)[4]+5)

#
system.time(BatMetrics(wave=smpl, info=PI, typeOfAnalysis = "B", noiseHet = 150,
                       myWL=512, amp.chk=5000))

#

### Analyse en série
file.choose()
setwd("/media/jeff/Médiathèque/Audio/Chiroptera/00_Backup_Recordings/0_Analyse_OK/20160712_LIFE_Rochefort")
list.files(pattern = "wav",ignore.case=T)
AllFiles <- list.files(pattern = "wav",ignore.case=T)[c(14,15)]
#filename <- list.files(pattern = "wav",ignore.case=T)[14]
PI <- "Heterodyne_BatMetrics_1"
sapply(AllFiles, function(x){
  filename <<- x
  smpl <- readWave(filename) # tuneR::readWave()
  BatMetrics(wave=smpl, info=PI, typeOfAnalysis = "H", noiseHet = 20,
             myWL=512, amp.chk=5000)
})

###################################
system.time(BatMetrics(wave=smpl, info="TrainMYODAU_faibles&forts", typeOfAnalysis = "T", myWL=256))


system.time(BatMetrics(wave=smpl, info=PI, typeOfAnalysis = "H"))
#BatMetrics(wave=smpl, info="TestHet", typeOfAnalysis = "H", noiseHet = noise)
#BatMetrics(wave=smpl, info="test", typeOfAnalysis = "B", myWL=256)

#
##
#####











#####
##
#
BatMetrics <- function(wave, info, typeOfAnalysis = c("H","T","B"),
                       myWL=512, amp.chk=1300, thSig=700, noiseHet=noise){
  ### PACKAGES ##########
  require('seewave')
  require('xlsx')
  require('tuneR')
  require('ggplot2')
  ### Points de repères ...
  Dur.Sig <- c(0.0008,0.020) # Pas utilisé
  Dur.Int <- c(0.010,0.600) # 600 ms = longest interval between calls
  ### HETERODYNE ####
  if(typeOfAnalysis == "H" | typeOfAnalysis == "B"){
    smp.lch <- channel(wave, "left")
    main.rg <- range(smp.lch@left)
    if((main.rg[2]-main.rg[1]) > 2000){ # Value of 2000 chosen empirically!! Discar empty samples
      png(paste0(substr(filename,1,nchar(filename)-4),"_GlobHet.png"),
          width = 600, height = 400, res=80)
      timer.smpl.peaks <- timer.abs(smp.lch, threshold=noiseHet,
                                    msmooth=c(220.5, 90),plotoutline = F) # Résolution à 0.005 seconde
      dev.off()
      #length(timer.smpl.peaks$s) # number of peaks detected
      #timer.smpl.peaks # detail
      
      if(timer.smpl.peaks$first == "pause"){
        i <- timer.smpl.peaks$p
        #        d <- timer.smpl.peaks$s
      } else {
        i <- timer.smpl.peaks$p[-1]
        timer.smpl.peaks$s.start <-  timer.smpl.peaks$s.start[-1]
        ##### !!!!! adapter pour d aussi cf. buzz
      }
      lg <- length(i)
      ip1 <- c(i[2:lg],0)
      ip2 <- c(i[3:lg],0,0)
      ip3 <- c(i[4:lg],0,0,0)
      ip4 <- c(i[5:lg],0,0,0,0)
      im1 <- c(0,i[1:lg-1])
      im2 <- c(0,0,i[2:lg-2])
      im3 <- c(0,0,0,i[3:lg-3])
      im4 <- c(0,0,0,0,i[4:lg-4])
      
      #      logi.ser <- (i < Dur.Int[2] & ip1 < Dur.Int[2]&
      #                     ip2 < Dur.Int[2] &
      #                     ip3 < Dur.Int[2] &
      #                     ip4 < Dur.Int[2] ) |
      #        (i < Dur.Int[2] & im1 < Dur.Int[2] &
      #           im2 < Dur.Int[2] &
      #           im3 < Dur.Int[2] &
      #           im4 < Dur.Int[2] ) ## Detect 5 following gaps of < 0.6 sec
      
      logi.ser <- ((i < Dur.Int[2] & i > Dur.Int[1]) &
                     (ip1 < Dur.Int[2] & ip1 > Dur.Int[1]) &
                     (ip2 < Dur.Int[2] & ip2 > Dur.Int[1]) &
                     (ip3 < Dur.Int[2] & ip3 > Dur.Int[1]) &
                     (ip4 < Dur.Int[2] & ip4 > Dur.Int[1])) |
        ((i < Dur.Int[2] & i > Dur.Int[1]) &
           (im1 < Dur.Int[2] & im1 > Dur.Int[1]) &
           (im2 < Dur.Int[2] & im2 > Dur.Int[1]) &
           (im3 < Dur.Int[2] & im3 > Dur.Int[1]) &
           (im4 < Dur.Int[2] & im4 > Dur.Int[1])) ## Detect 5 following gaps of < 0.6 sec
      logi.buz <- (i < Dur.Int[1] & ip1 < Dur.Int[1] & ip2 < Dur.Int[1] &
                     ip3 < Dur.Int[1] & ip4 < Dur.Int[1]) ## Detect 5 following gaps of < 0.01 sec
      
      PotBuz <- round(timer.smpl.peaks$s.start[logi.buz],2) ## Collect signals belonging to potential buzz
      signOK <- timer.smpl.peaks$s.start[logi.ser]
      #signOK # signals recognized as true
      brk <- which(signOK[2:length(signOK)] - signOK[1:length(signOK)-1] > 5)
      ## output the values of 5 full seconds as in the french method
      DF <- if(length(brk) == 0){
        rg <- range(signOK)[2] - range(signOK)[1]
        div5 <- rg %/% 5
        res0 <- ifelse(div5 == 0, 1,div5)
        c(Folder=fol.name,
          Full5sec=res0,
          TimeStart=round(signOK[1],3))} else {
            brk <- c(0,brk,length(signOK))
            sapply(1:(length(brk)-1), function(k){
              fst <- brk[k]+1; lst <- brk[k+1]
              rg <- range(signOK[fst:lst])[2] - range(signOK[fst:lst])[1]
              div5 <- rg %/% 5
              res <- ifelse(div5 == 0, 1,div5)
              c(Folder=fol.name,
                Full5sec=as.numeric(res),
                TimeStart=round(signOK[fst],3),
                TS.MinSec=paste0(signOK[fst] %/% 60,"'",round(signOK[fst] %% 60,1),"''"))
            })
          }
      DF.Het <- data.frame(t(DF))
      DF.Het$Full5sec <- as.numeric(as.character(DF.Het$Full5sec))
      DF.Het$TimeStart <- as.numeric(as.character(DF.Het$TimeStart))
    } else {DF.Het <- data.frame(Folder=fol.name, Full5sec= "Pas de contact hétérodyne")}
  }
  ### TIME EXPANSION ####
  if(typeOfAnalysis == "T" | typeOfAnalysis == "B"){
    smp.rch <- channel(wave, "right")
    ##### Detect chunks of time expansion
    png(paste0(substr(filename,1,nchar(filename)-4),"_GlobTE.png"),
        width = 600, height = 400, res=80)
    #tr <- timer(smp.rch, threshold=7, msmooth = c(4410,0), power=0.5) #resolution:1sec
    tr <- timer.abs(smpl@right, f=44100, threshold=10, msmooth = c(8820,50)) #resolution: 2 sec
    dev.off()
    # Correct timer to not include gaps of less tha 1.7 sec
    logiGapTE <- c(TRUE, tr$p[-c(1, length(tr$p))] > 1.7 ,TRUE)
    S_Start <- tr$s.start[logiGapTE[-length(logiGapTE)]]
    S_End <- tr$s.end[logiGapTE[-1]]
    tr$s.start <- S_Start
    tr$s.end <- S_End
    #if(tr$first == "signal"){
    #  s <- tr$s
    #  p <- tr$p
    #  s.start <- tr$s.start
    #  s.end <- tr$s.end
    #  tr$s <- p
    #  tr$p <- s[-1]
    #  tr$s.start <- s.end
    #  tr$s.end <- s.start[-1]
    #}
    Nchunk <- length(tr$s.start)
    names.chk <- paste0("Chunk_",1:Nchunk)
    ### make equal length of start and end vectors by adding -if unequal- total length
    if(Nchunk > length(tr$s.end)){
      tr$s.end[Nchunk] <- length(smp.rch) / smp.rch@samp.rate}
    chnk <- list()
    tr.chk <- list()
    logi.chk <- c()
    rangeTE <- c()
    thSigAuto <- numeric()
    ### Loop detection of signals
    for(j in 1:Nchunk){
      chnk[[j]] <- cutw(smp.rch,
                        from= tr$s.start[j], to=tr$s.end[j])
      mainTE.rg <- range(chnk[[j]])
      rangeTE[j] <- mainTE.rg[2] - mainTE.rg[1]
      ## Collect logical for actually analysed chunks
      logi.chk[j] <- rangeTE[j] > amp.chk
      if(logi.chk[j]){
        #        tr.chk[[j]] <- timer(chnk[[j]], f=44100, threshold=thSig, msmooth = c(44.1,0), power=1.2, plot=T)
        ## FIRST step detect the background noise. It calculates enveloppe and later the timer will
        ## reuse the same enveloppe. For the test phase enveloppe is computed two times, later I'll modify
        ## timer.abs() in order to RE-USE the enveloppe and speed up processing
        envel <- env(chnk[[j]], f=44100, envt="abs", msmooth = c(44.1,90), plot=F)
        name.diff_1 <- names(which(abs(diff(table(cut(envel,1000)))) < 1000)[1])
        thSigAuto[j] <- as.numeric(substring(name.diff_1,2,regexpr(".",name.diff_1, fixed=T)-1)) + 1 # overrides the thSig provided as parameter in the function. If it works well, thSig will be deprecated
        tr.chk[[j]] <- timer.abs(chnk[[j]], f=44100, threshold=thSigAuto[j], msmooth = c(44.1,90), plot=T)
        ### make equal length of start and end vectors by adding -if unequal- total length
        if(length(tr.chk[[j]]$s.start) > length(tr.chk[[j]]$s.end)){
          tr.chk[[j]]$s.end[length(tr.chk[[j]]$s.start)] <- (length(chnk[[j]]) / 44100)-0.05}
      }
    }
    names(chnk[logi.chk]) <- names(tr.chk[logi.chk]) <- names.chk[logi.chk]
    #paste(paste(names.chk,rangeTE,sep=": "), collapse=" ; ")
    #### Measurments of all chunks
    lis <- lapply((1:Nchunk)[logi.chk], function(cnb){
      if(class(tr.chk[[cnb]]) == "list"){
        ############## Edit the result of timer() to solve signals cut into 2 because of narrowing of intensity
        if(length(tr.chk[[cnb]]$p) > length(tr.chk$s)){tr.chk[[cnb]]$p <- tr.chk[[cnb]]$p[1:length(tr.chk[[cnb]]$s)]}
        logiNotSplitSignal <- !(tr.chk[[cnb]]$p < 0.003 &
                                  c(0,head(tr.chk[[cnb]]$s,-1)) > 0.01 &
                                  c(tr.chk[[cnb]]$s[-1],0) > 0.01)
        tr.chk[[cnb]]$p <- tr.chk[[cnb]]$p[logiNotSplitSignal]
        tr.chk[[cnb]]$s.start <- tr.chk[[cnb]]$s.start[logiNotSplitSignal]
        tr.chk[[cnb]]$s.end <- tr.chk[[cnb]]$s.end[c(logiNotSplitSignal[-1],T)]
        tr.chk[[cnb]]$s <- tr.chk[[cnb]]$s.end - tr.chk[[cnb]]$s.start
        ##############
        logi.len <- tr.chk[[cnb]]$s > 0.010 # ~min threshold 0.011 for signals (!! may ignore some short MYOBRA !!)~
        Nsig <- length(tr.chk[[cnb]]$s[logi.len])
        st <- tr.chk[[cnb]]$s.start[logi.len]
        nd <- tr.chk[[cnb]]$s.end[logi.len]
        lgt <- tr.chk[[cnb]]$s[logi.len]
        ## Fix problem when st = 0 !?
        logi.No0 <- st != 0
        st <- st[logi.No0]
        nd <- nd[logi.No0]
        lgt <- lgt[logi.No0]
        Nsig <- Nsig - sum(!logi.No0)
        ## Fix case when st > nd
        logiSt.nd.FALSE <- st < nd
        st <- st[logiSt.nd.FALSE]
        nd <- nd[logiSt.nd.FALSE]
        lgt <- lgt[logiSt.nd.FALSE]
        Nsig <- Nsig - sum(!logiSt.nd.FALSE)
        if(Nsig > 0){
          # Measurments inside of 1 chunk
          fme <- fi <- ft <- ft1 <- ft2 <- numeric()
          for(subchk in 1:Nsig){
            # measure FME
            sp.peak <- spec(chnk[[cnb]], f=44100, wl=myWL, plot=F,
                            from = st[subchk],to = nd[subchk])
            fme[subchk] <- sp.peak[which.max(sp.peak[,2]),1]
            # measure FI
            sp.peak <- spec(chnk[[cnb]], f=44100, wl=myWL, plot=F,
                            at = st[subchk])
            fi[subchk] <- sp.peak[which.max(sp.peak[,2]),1]
            # measure FT --> NOTE it checks if the value at t_end is smaller (cf. QFC)!!
            sp.peak.1 <- spec(chnk[[cnb]], f=44100, wl=myWL, plot=F,
                              at = nd[subchk])
            sp.peak.2 <- spec(chnk[[cnb]], f=44100, wl=myWL, plot=F,
                              at = nd[subchk]-0.01) ### still useful end-0.01 ??
            ft1[subchk] <- sp.peak.1[which.max(sp.peak.1[,2]),1]
            ft2[subchk] <- sp.peak.2[which.max(sp.peak.2[,2]),1]
            ft[subchk] <- min(c(ft1[subchk],ft2[subchk]))
          }
          DF <- data.frame(Chunk = names.chk[cnb],
                           FME = round(fme*10, 1),
                           FI = round(fi*10, 1),
                           FT = ifelse(fme > ft,round(ft*10, 1),round(fme*10, 1)), #take the smallest value between FME and FT (QFC!)
                           LB = ifelse(fme > ft,round((fi-ft)*10 ,1), round((fi-fme)*10 ,1)), #take the smallest value between FME and FT (QFC!)
                           Durée = round(lgt*100, 2),
                           StartInSeq = round(st+tr$s.start[cnb],3))
          #          if(Nsig > 1){itv <- 100 * (c(DF$StartInSeq[-1],0) - DF$StartInSeq)
          #          DF$Intervalle <- c(itv[-dim(DF)[1]], NA)} else {DF$Intervalle <- NA}
          DF$DiagQFC <- ifelse(DF$LB < 5,"QFC?","-")
          DF}} else {DF <- data.frame(Chunk = names.chk[cnb],
                                      FME = NA,
                                      FI = NA,
                                      FT = NA,
                                      LB = NA,
                                      Durée = NA,
                                      StartInSeq = NA,
                                      Intervalle = NA,
                                      DiagQFC = NA)}
    })
    #lis
    DF.TE <- do.call("rbind",lis)
    if(dim(DF.TE)[1] != 0){
      if(dim(DF.TE)[1] > 1){itv <- 100 * (c(DF.TE$StartInSeq[-1],0) - DF.TE$StartInSeq)
      DF.TE$Intervalle <- c(itv[-dim(DF.TE)[1]], NA)} else {DF.TE$Intervalle <- NA}
    }    
    discNul <- DF.TE$Durée < 5 & DF.TE$LB < 7 | # Discard echoes (e.g.on QFC) and voice
      DF.TE$Intervalle < 5 |
      DF.TE$Durée > 25 | #discard all Durée > 35 ms, extremes!!?? changed to 25 
      DF.TE$FME < 10 | DF.TE$FME > 90 | # discard extremes FME
      DF.TE$FI < 25 & DF.TE$Durée < 7 | # discard low FI which are not QFC
      DF.TE$FT < 15 & DF.TE$Durée < 9 | # discard low FT which are not QFC
      DF.TE$FT < 6 | # discard all FT lower than 6 kHz!
      DF.TE$FT > 50 | # discard high FT
      DF.TE$FME > DF.TE$FI # discard impossible values i.e. FME > FI
    DF.TE <- DF.TE[!discNul,]
    if(dim(DF.TE)[1] != 0){
      if(dim(DF.TE)[1] > 1){itv <- 100 * (c(DF.TE$StartInSeq[-1],0) - DF.TE$StartInSeq)
      DF.TE$Intervalle <- c(itv[-dim(DF.TE)[1]], NA)} else {DF.TE$Intervalle <- NA}
    }    
    #DF.TE
    #### GGPLOT GRAPHS ######
    ggplot(DF.TE, aes(x = Durée, y = FT, colour = Chunk)) +
      geom_point() + ggtitle(info) +
      geom_hline(yintercept = c(23,30), linetype=c(2,2), size=c(0.4,0.4)) +
      scale_y_continuous(limits=range(DF.TE$FT)) + theme_grey(12)
    ggsave(paste0(substr(filename,1,nchar(filename)-4),"_Plot_FT-Durée.png"),
           scale=1.5)
    ggplot(DF.TE, aes(x = LB, y = FME, colour = Chunk)) +
      geom_point() + ggtitle(info) +
      geom_vline(xintercept = 5, linetype=2, size=0.7) +
      scale_x_continuous(limits=range(DF.TE$LB)) + theme_grey(12)
    ggsave(paste0(substr(filename,1,nchar(filename)-4),"_Plot_FME-LB.png"),
           scale=1.5)
  }
  #### EXPORT XLSX ######
  wb <- createWorkbook(type="xlsx")
  # Styles
  cs_blue <- CellStyle(wb) + Font(wb, heightInPoints=10, color="blue")  # blue caract
  cs_head <- CellStyle(wb) + Font(wb, heightInPoints=10, isBold=TRUE) + Border() + Alignment(h="ALIGN_CENTER") # header
  cs_numTab <- CellStyle(wb) + Font(wb, heightInPoints = 10) + Alignment(h="ALIGN_CENTER")
  cs_small <- CellStyle(wb) + Font(wb, heightInPoints = 8) + Alignment(h="ALIGN_CENTER")
  SHEET_x <- createSheet(wb,paste0(info,"_WL",myWL))
  # Feuille 1
  if(typeOfAnalysis == "H" | typeOfAnalysis == "B"){
    addDataFrame(DF.Het,SHEET_x
                 ,col.names=TRUE, row.names=TRUE
                 ,startRow=1,startColumn=1
                 ,colnamesStyle=cs_head, 
                 rownamesStyle=cs_blue,
                 colStyle=list(`2`=cs_numTab, `3`=cs_numTab, `4`=cs_numTab, `5`=cs_numTab,
                               `6`=cs_numTab, `7`=cs_numTab, `8`=cs_numTab))
    if(exists("PotBuz")){addDataFrame(data.frame(PotBuz), SHEET_x, col.names = T, row.names = F,
                                      startColumn = 1, startRow = nrow(DF.Het) + ifelse(typeOfAnalysis == "H",20,30),
                                      colnamesStyle=cs_small,
                                      colStyle=list(`1`=cs_small, `2`=cs_small, `3`=cs_small))}
    autoSizeColumn(SHEET_x, 1:5)
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_GlobHet.png"), scale=0.5
               ,sheet=SHEET_x,startColumn = 1, startRow = nrow(DF.Het)+3)
  }
  if(typeOfAnalysis == "B"){
    addDataFrame(DF.TE,SHEET_x
                 ,col.names=TRUE, row.names=TRUE
                 ,startRow=1,startColumn=7
                 ,colnamesStyle=cs_head,
                 rownamesStyle=cs_blue,
                 colStyle=list(`2`=cs_numTab, `3`=cs_numTab, `4`=cs_numTab, `5`=cs_numTab,
                               `6`=cs_numTab, `7`=cs_numTab, `8`=cs_numTab))
    autoSizeColumn(SHEET_x, 6:17)
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_GlobTE.png"), scale=0.5
               ,sheet=SHEET_x,startColumn = 1, startRow = nrow(DF.Het)+16)
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_Plot_FT-Durée.png"), scale=0.7
               ,sheet=SHEET_x,startColumn = 17)
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_Plot_FME-LB.png"), scale=0.7
               ,sheet=SHEET_x,startColumn = 25)
    addDataFrame(data.frame(Chunk = names.chk,
                            Amplitude = rangeTE,
                            Accepted = logi.chk,
                            AutoThsld= thSigAuto,
                            Start=paste0(tr$s.start %/% 60,"'",round(tr$s.start %% 60,1),"''")),
                 SHEET_x, col.names = T, row.names = F,
                 startColumn = 2, startRow = nrow(DF.Het)+30,
                 colnamesStyle=cs_small,
                 colStyle=list(`1`=cs_small, `2`=cs_small, `3`=cs_small))
    if(exists("PotBuz")){addDataFrame(data.frame(PotBuz), SHEET_x, col.names = T, row.names = F,
                                      startColumn = 1, startRow = nrow(DF.Het)+30,
                                      colnamesStyle=cs_small,
                                      colStyle=list(`1`=cs_small, `2`=cs_small, `3`=cs_small))}
  }
  
  if(typeOfAnalysis == "T"){
    addDataFrame(DF.TE,SHEET_x
                 ,col.names=TRUE, row.names=TRUE
                 ,startRow=1,startColumn=1
                 ,colnamesStyle=cs_head,
                 rownamesStyle=cs_blue,
                 colStyle=list(`2`=cs_numTab, `3`=cs_numTab, `4`=cs_numTab, `5`=cs_numTab,
                               `6`=cs_numTab, `7`=cs_numTab, `8`=cs_numTab))
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_GlobTE.png"), scale=0.7
               ,sheet=SHEET_x,startColumn = 1, startRow = nrow(DF.TE)+3)
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_Plot_FT-Durée.png"), scale=0.8
               ,sheet=SHEET_x,startColumn = 13)
    addPicture(file=paste0(substr(filename,1,nchar(filename)-4),"_Plot_FME-LB.png"), scale=0.8
               ,sheet=SHEET_x,startColumn = 19)
    addDataFrame(data.frame(Chunk = names.chk,
                            Amplitude = rangeTE,
                            Accepted = logi.chk,
                            Start=paste0(tr$s.start %/% 60,"'",round(tr$s.start %% 60,1),"''")),
                 SHEET_x, col.names = T, row.names = F,
                 startColumn = 2, startRow = nrow(DF.TE)+21,
                 colnamesStyle=cs_small,
                 colStyle=list(`1`=cs_small, `2`=cs_small, `3`=cs_small))
    autoSizeColumn(SHEET_x, 1:10)
  }
  filenameNoExt <-  substring(filename,1,nchar(filename)-4)
  saveWorkbook(wb, paste0(fol.name,"_",filenameNoExt,"_",info,"_WL",myWL,".xlsx"))
  
}

#### The function timer{seewave} of Jérôme Sueur is adapted for my needs to control an absolute
#### scale allowing to define a common noise background and thresholds for signals detection.
#### I also modified the 'plot' part of the function.

timer.abs <- function (wave, f, threshold = 400, dmin = NULL, envt = "abs", 
                       power = 1, msmooth = NULL, ksmooth = NULL, ssmooth = NULL, 
                       tlim = NULL, plot = TRUE, plotthreshold = TRUE,
                       plotoutline=TRUE, col = "black", 
                       colval = "red", xlab = "Time (s)", ylab = "Amplitude", ...) 
{
  input <- inputw(wave = wave, f = f)
  wave <- input$w
  f <- input$f
  rm(input)
  n <- length(wave)
  thres <- threshold
  if (power == 0) 
    stop("'power' cannot equal to 0")
  if (!is.null(msmooth) && !is.null(ksmooth)) 
    stop("'msmooth' and 'ksmooth' cannot be used together")
  if (!is.null(msmooth) && !is.null(ssmooth)) 
    stop("'msmooth' and 'ssmooth' cannot be used together")
  if (!is.null(ksmooth) && !is.null(ssmooth)) 
    stop("'ksmooth' and 'ssmooth' cannot be used together")
  if (!is.null(dmin)) {
    if (length(dmin) != 1) 
      stop("'dmin' should be a numeric vector of length 1")
    if (dmin <= 0) 
      stop("'dmin' cannot be negative or equal to 0")
    if (dmin >= n/f) 
      stop("'dmin' cannot equal or be higher than the wave duration")
  }
  if (!is.null(tlim)) {
    wave <- cutw(wave, f = f, from = tlim[1], to = tlim[2])
    n <- length(wave)
  }
  wave1 <- env(wave = wave, f = f, msmooth = msmooth, ksmooth = ksmooth, 
               ssmooth = ssmooth, envt = envt, norm = FALSE, plot = FALSE)
  n1 <- length(wave1)
  f1 <- f * (n1/n)
  if (power != 1) 
    wave1 <- wave1^power
  wave2 <- ifelse(wave1 <= thres, yes = 1, no = 2)
  n2 <- length(wave2)
  wave4 <- apply(as.matrix(1:(n2 - 1)), 1, function(x) wave2[x] + 
                   wave2[x + 1])
  n4 <- length(wave4)
  wave4[c(1, n4)] <- 3
  wave5 <- which(wave4 == 3)
  if (!is.null(dmin)) {
    event.dur <- diff(wave5)
    event.idx <- which(event.dur < dmin * f1)
    if (length(event.idx) != 0) {
      for (i in event.idx) {
        wave4[(wave5[i]):(wave5[i] + event.dur[i])] <- 2
      }
      wave4[which(abs(diff(wave4)) == 2)] <- 3
      wave4[c(1, n4)] <- 3
    }
    wave5 <- which(wave4 == 3)
    if (length(wave5) == 2) {
      stop("'dmin' was set to a too high value, there are no signal longer than 'dmin'")
    }
  }
  wave5[-1] <- wave5[-1] + 1
  f4 <- f * (n4/n)
  wave4 <- ts(wave4, start = 0, end = n4/f4, frequency = f4)
  positions <- time(wave4)[wave5]
  durations <- diff(positions)
  npos <- length(positions)
  if (wave2[1] == 1) {
    first <- "pause"
    pause <- durations[seq(1, npos - 1, by = 2)]
    signal <- durations[seq(2, npos - 1, by = 2)]
    start.signal <- positions[seq(2, npos - 1, by = 2)]
    end.signal <- positions[seq(3, npos - 1, by = 2)]
  }
  else {
    first <- "signal"
    pause <- durations[seq(2, npos - 1, by = 2)]
    signal <- durations[seq(1, npos - 1, by = 2)]
    start.signal <- positions[seq(1, npos - 1, by = 2)]
    end.signal <- positions[seq(2, npos - 1, by = 2)]
  }
  ratio <- sum(signal)/sum(pause)
  timer <- list(s = signal, p = pause, r = ratio, s.start = start.signal, 
                s.end = end.signal, first = first)
  if (plot) {
    plot(x = seq(0, n1/f1, length.out = n1), y = wave1, xlab = xlab, 
         ylab = ylab, ylim = c(0, 1.1*max(wave1)), col = col, #, yaxt = "y"
         type = "l", xaxs = "i", ...)
    if (plotthreshold) {
      abline(h = thres, col = colval, lty = 2)
      mtext(as.character(threshold), side = 2, 
            line = 0.5, at = thres, las = 1, col = colval, 
            cex = 0.8)
    }
    if(plotoutline){
      outline <- wave4 - 3
      outline[outline < 0] <- 0
      outline[1] <- NA
      lines(x = seq(0, n1/f1, length.out = n1), y = c(outline)*max(wave1), 
            col = colval)
      wave8 <- numeric(npos - 1)
      for (i in 2:npos) {
        wave8[i] <- ((wave5[i] - wave5[i - 1])/2) + wave5[i - 
                                                            1]
      }
      if (wave2[1] == 1) {
        wave8.1 <- wave8[seq(2, npos, by = 2)]/f1
        wave8.2 <- wave8[seq(3, npos, by = 2)]/f1
      }
      else {
        wave8.2 <- wave8[seq(2, npos, by = 2)]/f1
        wave8.1 <- wave8[seq(3, npos, by = 2)]/f1
      }
      ypl <- as.character(round(pause, 2))
      ysl <- as.character(round(signal, 2))
      text(x = wave8.1, y = 0.075*max(wave1), ypl, col = colval, cex = 0.8)
      text(x = wave8.2, y = 1.075*max(wave1), ysl, col = colval, cex = 0.8)
    }
    invisible(timer)
  }
  else {
    return(timer)
  }
}

