#' ---
#' title: 'Patients with presumed tuberculosis in sub-Saharan Africa that are not diagnosed with tuberculosis: a systematic review and meta-analysis (statistical appendix)'
#' author:
#' - S Jayasooriya,
#' - F Dimambro-Denson,
#' - C Beecroft,
#' - J Balen,
#' - B Awokola,
#' - C Mitchell,
#' - B Kampmann,
#' - F Campbell,
#' - PJ Dodd,
#' - K Mortimer
#' date: "August, 2021"
#' output:
#'    pdf_document:
#'      toc: true
#'      highlight: kate
#' ---

#' # Pre-amble
#'
#' This document is generated from an R script in literate programming fashion. All R code is quoted in this document, together with output (preceeded by '##') and figures. The article forest plot is saved to the `output` folder but not included in the document since it is too cramped. The script and data are publicly available on GitHub at  https://github.com/petedodd/NotTB and once the repository is downloaded, it should be possible to generate this document using R with the command
#' `rmarkdown::render('NotTBmeta.R')` within R, or from a unix-like command line with `R -q -e "rmarkdown::render(\"NotTBmeta.R\",output_dir=\"./output\")"`. Alternatively, the R script can be run in whole or part as a conventional R script.
#'
#' ## Dependencies
#'
#' To compile this document, the rmarkdown & knitr packages must be installed.
#' The other R packages required to run this analysis should be installed if necessary, and loaded, with:
#'
pkgs.needed <- c("ggplot2","scales","cowplot","ggpubr", #graphs
                 "data.table","here",                   #data mgt
                 "metafor")                             #metaanalysis
install.packages(setdiff(pkgs.needed, rownames(installed.packages())))
suppressMessages(
    devnull <- lapply(pkgs.needed, require, character.only = TRUE) #load for use
)

#'
#' This analysis was run using:
#'
sI <- sessionInfo()
dI <- data.frame(
    item=c('R version','platform','OS','metafor version'),
    version=c(
        sI$R.version$version.string, #R version
        sI$platform,                 #platofm
        sI$running,                  #OS
        sI$otherPkgs$metafor$Version #metafor version
    )
)
knitr::kable(dI)


#'
#' # Main analyses
#'
#' ## Approach
#' We use a generalized linear mixed effects (GLMM) approach to meta-analysis assuming a binomial response and logit link^[[Stijnen T, Hamza TH, Ozdemir P. Random effects meta-analysis of event outcome in the framework of the generalized linear mixed model with applications in sparse data. ](http://dx.doi.org/10.1002/sim.4040)]. This means we assume
#'
#' $$k_i\sim \mathrm{Binomial}(N_i,p_i)$$
#' $$\mathrm{logit}(p_i) = \mu+\varepsilon_i $$
#' $$ \varepsilon_i \sim \mathcal{N}(0,\tau^2)$$
#'
#' where $i=1,\dots,S$ indexes the numbers of studies.
#'
#' Use of arcsine or double arcsine transformations has been criticized in this context, with the GLMM.^[[Schwarzer G, Chemaitelly H, Abu-Raddad LJ, Rücker G. Seriously misleading results using inverse of Freeman-Tukey double arcsine transformation in meta-analysis of single proportions](https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1348)]
#' 
#'

#' Read in the data and ensure that factors behave as intended:
DD <- fread(file=here('SRMAdata.csv'))
DD[,lab:=factor(lab,levels=rev(DD[order(bac)]$lab),ordered = TRUE)]

#' Create exact binomial confidence intervals:
ciz <- function(x,y){
    x <- as.integer(x); y <- as.integer(y)
    list(binom.test(x,y)$conf.int[1],binom.test(x,y)$conf.int[2])
}
DD[,`NotTB Proportion`:=NnotTB/N]
for(i in 1:nrow(DD)){ DD[i,c('lo','hi'):=ciz(NnotTB,N)]; }
DD[,SE:=(hi-lo)/3.92]

#'
#' ## Meta-analyses
#' 
#' Meta-analysis for passively found TB patients with bacteriologically unconfirmed TB included:
maPU <- rma.glmm(measure = "PLO", #  binomial w/ logit link
            xi = NnotTB,     # numerator
            ni = N,          # denominator
            data = DD[mode=='Passive' &
                      clinical=='(Unconfirmed TB included)'],
            slab = Author)      # what to use as labels on graphs
summary(maPU)
forest(maPU,transf = transf.ilogit,refline=NA)

#' Meta-analysis for passively found TB patients with bacteriologically unconfirmed TB excluded:
maPN <- rma.glmm(measure = "PLO", #  binomial w/ logit link
            xi = NnotTB,     # numerator
            ni = N,          # denominator
            data = DD[mode=='Passive' &
                      clinical=='(No unconfirmed TB)'],
            slab = Author)      # what to use as labels on graphs
summary(maPN)
forest(maPN,transf = transf.ilogit,refline=NA)

#' Meta-analysis for actively found TB patients:
maA <- rma.glmm(measure = "PLO", #  binomial w/ logit link
            xi = NnotTB,     # numerator
            ni = N,          # denominator
            data = DD[mode=='Active'],
            slab = Author)      # what to use as labels on graphs
summary(maA)
forest(maA,transf = transf.ilogit,refline=NA)

#' Make predictions for plot data:
map <- predict(maA,transf = transf.ilogit)
mup <- predict(maPU,transf = transf.ilogit)
mnp <- predict(maPN,transf = transf.ilogit)


#'
#' ## Creation of combined forest plot
#' 
#' Summary data for combined forest plot:
f1 <- function(x)format(round(x,1),nsmall=1)
cnz <- c('(Unconfirmed TB included)',
         '(No unconfirmed TB)',
         '(No unconfirmed TB)')
predz <- data.table(mode=c('Passive','Passive','Active'),
                    clinical=cnz,
                    `NotTB Proportion` = c(mup$pred,mnp$pred,map$pred),
                    lo = c(mup$ci.lb,mnp$ci.lb,map$ci.lb),
                    hi = c(mup$ci.ub,mnp$ci.ub,map$ci.ub),
                    lab=paste0('SUMMARY (',expression(I^2),'=',
                               f1(c(maA$I2,maPN$I2,maPU$I2)),'%)')
                    )
predz[,SE:=(hi-lo)/3.92]
predz[,qty:='summary']
predz[,bac:=0]
predz[,mid:=`NotTB Proportion`]
predz[,CI:=paste0(f1(1e2*mid),' (',f1(1e2*lo),' - ',f1(1e2*hi),')')]
predz[,wt:='100.0%']
predz[,w:=1]

#' Appending plot data to inputs:
DD[,qty:='study']
DD[,mid:=`NotTB Proportion`]
DD[,CI:=paste0(f1(1e2*mid),' (',f1(1e2*lo),' - ',f1(1e2*hi),')')]
DD[,wt:=1/SE^2]
DD[,wtt:=sum(wt),by=.(mode,clinical)]
DD[,wt:=1e2*wt/wtt]
DD[,wt:=paste0(f1(wt),'%')]
DD[,w:=0]

#' Combined plot data:
B <- rbind(
    DD[,.(lab,`NotTB Proportion`,lo,hi,SE,mode,clinical,
          qty,bac,CI,wt,w)],
    predz[,.(lab,`NotTB Proportion`,lo,hi,SE,mode,clinical,
             qty,bac,CI,wt,w)]
)
lbz <- as.character(B[order(bac)]$lab)
lbz2 <- c(lbz[1:3],rev(lbz[-c(1:3)]))
B[,lab:=factor(lab,levels=lbz2,ordered = TRUE)]
B[,clinical.g:='Clinically diagnosed tuberculosis included']
B[clinical=='(No unconfirmed TB)',
  clinical.g:='No clinically diagnosed tuberculosis included']
B[mode=='Active',clinical.g:='']
B[,mode:=factor(mode,levels=c('Passive','Active'),ordered = TRUE)]
B[,clinical.g:=factor(clinical.g,levels=unique(clinical.g))]
labdat <- B[1]
labdat[,txt:=' weight (%)']

#' Create publication forest plot figure:
SA <- ggplot(B,aes(lab,y=`NotTB Proportion`,
                   ymin=lo,
                   ymax=hi,
                   col=qty)) +
    geom_point(aes(size=1/SE^2,shape=qty)) +
    geom_errorbar(aes(width=w/2)) +
    scale_y_continuous(label=percent,limits = c(0,NA))+
    scale_color_manual(values=c('study'="black",'summary'="blue"))+
    scale_shape_manual(values=c('study'=22,'summary'=23))+
    xlab('') +
    ylab('Proportion of patients with presumptive tuberculosis not diagnosed as tuberculosis')+
    facet_grid(mode + clinical.g ~ .,
               scales = 'free',space='free',
               switch='x'
               )+
    coord_flip() +
    guides(size='none',color='none',shape='none')+
    theme_classic() +
    theme(panel.spacing = unit(2, "lines"), #or 3
          strip.background = element_blank(),
          strip.placement = "outside") +
    geom_text(aes(x=lab,y=1.2,label=CI,hjust='right')) +
    geom_text(aes(x=lab,y=0.0,label=wt))+
    geom_text(data=labdat,aes(x=9.5,y=0,label=txt))+
    ggpubr::grids()

ggsave(SA,file=here('output/ForestPlot.pdf'),h=13,w=12)
ggsave(SA,file=here('output/ForestPlot.eps'),h=13,w=12)


#'
#'
#' # Meta-regressions
#'
#' In this section we consider various potential sources of heterogeneity through scatter plots and meta-regression.
#'
#' ## TB prevalence
#'
#' The burden of TB in a population might reasonably be expected to influence the proportion of presumptive TB that is not TB.
#'
DD[,tb:=`WHO TB estimate (per 100 000 year of study)`]

ggplot(DD,aes(tb,`NotTB Proportion`,
              size=N,col=mode,shape=clinical))+
    scale_x_continuous(label=comma,limits=c(0,NA))+
    scale_y_continuous(label=percent,limits=c(0,1))+
    geom_point()+
    xlab('WHO estimate of TB prevalence per 100,000 for country-year')+
    ylab('Proportion not TB in study')+
    ggtitle('Influence of population TB burden')

#' We can formally investigating the influence of TB burden in explaining heterogeneity with a meta-regression:
tbmr <-  rma.glmm(measure = "PLO",  #binomial w/ logit link
              xi = NnotTB,     # numerator
              ni = N,          # denominator
              data = DD,        # what data to use
              mods = ~mode*clinical + tb)
summary(tbmr)


#'
#' ## HIV prevalence
#'
#' Population HIV prevalence may plausibly influence the proportion of presumptives not diagnosed with TB both by inluencing TB burden, but also by changing the typical clinical characteristics of TB and most importantly, the burden of other illness that could be designated presumptive TB.
#'
#' 
#'
ggplot(DD,aes(hiv/1e2,`NotTB Proportion`,
              size=N,col=mode,shape=clinical))+
    scale_x_continuous(label=percent,limits=c(0,0.13))+
    scale_y_continuous(label=percent,limits=c(0,1))+
    geom_point()+
    xlab('UNAIDS estimate of HIV prevalence 15-49 for country-year')+
    ylab('Proportion not TB in study')+
    ggtitle('Influence of population HIV prevalence')


#' We can formally investigating the influence of HIV in explaining heterogeneity with a meta-regression:
hivmr <-  rma.glmm(measure = "PLO",  #binomial w/ logit link
               xi = NnotTB,     # numerator
               ni = N,          # denominator
               data = DD,        # what data to use
               mods = ~mode*clinical + hiv)
summary(hivmr)

#'
#' ## Calendar time
#'
#' To explore whether there has been any change over time, we consider calendar year
#'
#'
ggplot(DD,aes(Year,`NotTB Proportion`,
              size=N,col=mode,shape=clinical))+
    scale_y_continuous(label=percent,limits=c(0,1))+
    geom_point()+
    xlab('Study year')+
    ylab('Proportion not TB in study')+
    ggtitle('Influence of calendar year')


#' We can formally investigating the influence of year in explaining heterogeneity with a meta-regression:
yearmr <-  rma.glmm(measure = "PLO",  #binomial w/ logit link
               xi = NnotTB,     # numerator
               ni = N,          # denominator
               data = DD,        # what data to use
               mods = ~mode*clinical + Year)
summary(yearmr)


#'
#' # Sensitivity analyses
#'
#' ## Dorman et al. by country only
#'
#' In the main analysis, we considered the different sites in the 2018 study by Dorman et al to be separate data. This included considering the two sites in South Africa - Cape Town and Johannesburg - as different, which was motivated by the very distinct TB epidemiology in the Western Cape. Here we investigate the impact of aggregating the two South African sites in Dorman et al on the meta-analysis for studies with passive case finding excluding clinically diagnosed TB.
#'
#' Restrict to relevant data & aggregate over Dorman in South Africa:
tmp <- DD[mode=='Passive' & clinical=='(No unconfirmed TB)']
tmp[,Country.Simple:=gsub(" \\-.+$","",Country)]                    #remove cities
tmp[,authorcountry:=paste(gsub("^([A-Za-z]+).*","\\1",Author),Country.Simple,sep = ", ")] #new label
tmp <- tmp[,.(NnotTB=sum(NnotTB),N=sum(N)),by=authorcountry]
knitr::kable(tmp) #check

#' Rerun this meta-analysis with the new data:
maPNsa <- rma.glmm(measure = "PLO", # binomial w/ logit link
            xi = NnotTB,     # numerator
            ni = N,          # denominator
            data =tmp,       # new data
            slab = authorcountry)      # what to use as labels on graphs
summary(maPNsa)
forest(maPNsa,transf = transf.ilogit,refline=NA)

#'
#'This is very similar to the main analysis above.
#'
#' ## Regional groupings
#'
#' Here we investigate whether country can explain some heterogeneity. Since when countries have occur only once, it is not possible to identify a country coefficient, we these countries into an "Other" category.
#'
DD[,Country.Group:=gsub(" \\-.+$","",Country)]                    #remove cities
DD[!Country.Group %in% c("South Africa","Ethiopia","Nigeria"),Country.Group:="Other"] #group
DD[,Country.Group:=factor(Country.Group,levels=unique(Country.Group))] #make factor

#' Plot this data:
ggplot(DD,aes(Country.Group,`NotTB Proportion`,
              size=N,col=mode,shape=clinical))+
    scale_y_continuous(label=percent,limits=c(0,1))+
    geom_point()+
    xlab('Country or country-group')+
    ylab('Proportion not TB in study')+
    ggtitle('Influence of region')


#' Perform meta-regression on country-group:
cgmr <-  rma.glmm(measure = "PLO",  #binomial w/ logit link
               xi = NnotTB,     # numerator
               ni = N,          # denominator
               data = DD,        # what data to use
               mods = ~mode*clinical + Country.Group)
summary(cgmr)

