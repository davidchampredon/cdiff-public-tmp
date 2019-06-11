args = commandArgs(trailingOnly = TRUE)
scen = args[1]

# DEBUG: scen=2

if(scen==1) rdataname <- 'simul-mc-frailty-40-prm_vax_p_0.RData'  #''
if(scen==2) rdataname <- 'simul-mc-frailty-80-prm_vax_p_0.RData' # ''

if(scen==1) rdataname.simple <- 'simul-mc-frailty-40-prm_vax_p_0-simplified.RData'  #'simul-mc-frailty-40-prm_vax_p_0.RData'
if(scen==2) rdataname.simple <- 'simul-mc-frailty-80-prm_vax_p_0-simplified.RData' # 'simul-mc-frailty-80-prm_vax_p_0.RData'

load(rdataname)  ;   
load(rdataname.simple)  ; df.s=simple.df ; rm(simple.df)

# a <- df[df$time>1825,]
library(dplyr);library(tidyr) #library(tidyverse)

print(rdataname)
nmc = length(unique(df.s$mc))
print(paste('Total MC iter: ',nmc))

print("  Susceptibility index of all simulations:")
summary(df$sIdx)

df.after = df[df$time>1825, ]
print(paste('proportion vax: ',mean(df.after$isVax)))


vaxeff_empirical <- function(i, df.s) {
    
    if(i>0) df.s = df.s[df.s$mc==i,]
    
    # Only colonized:
    a = df.s[df.s$wasColonized==TRUE,]
    
    # Among the vaccinated:
    n.col.vax   <- nrow(a[a$isVax==TRUE, ])
    n.col.s.vax <- nrow(a[a$totDurSymptoms>0 & a$isVax==TRUE, ])
    p.s.vax = n.col.s.vax / n.col.vax
    
    # Among not vaccinated:
    n.col.novax   <- nrow(a[a$isVax==FALSE, ])
    n.col.s.novax <- nrow(a[a$totDurSymptoms>0 & a$isVax==FALSE, ])
    p.s.novax = n.col.s.novax / n.col.novax

    # Vaccine efficacy    
    ve = 1 - p.s.vax/p.s.novax
    return(ve)
}

vv = sapply(X = 1:nmc, FUN = vaxeff_empirical, df.s=df.s)

print(paste("Empirical vax eff for",rdataname.simple,':'))
summary(vv)
pdf(paste0('plot-ve-emp-',scen,'.pdf'))
hist(vv, main=paste('Empirical vax eff\n',rdataname.simple), col='grey')
dev.off()

