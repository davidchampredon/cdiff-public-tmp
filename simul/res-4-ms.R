###
###  Results for manuscript
###

library(tidyr)
library(dplyr)

# ---- incidence -----

dat.inc <- read.csv('figure-1-data.csv') 

df.inc <- dat.inc %>%
    select(year.after.vax, m, scen) %>%
    spread(scen, m)
df.inc

inc.40 = df.inc$`simul-mc-frailty-40-prm_vax_p_0.RData`
inc.80 = df.inc$`simul-mc-frailty-80-prm_vax_p_0.RData`

chg.40 <- inc.40[3]/inc.40[1]-1
chg.80 <- inc.80[3]/inc.80[1]-1

print(df.inc)
print(paste('Reduction of annual incidence, scen_p_40 =', -chg.40))
print(paste('Reduction of annual incidence, scen_p_80 =', -chg.80))


# ---- isolations ----

dat.iso.40 <- read.csv('out-compare-iso-frailty-40-prm_vax_p_0.csv')
dat.iso.80 <- read.csv('out-compare-iso-frailty-80-prm_vax_p_0.csv')


print("Isolation days p_40")
dat.iso.40 %>%
    select(contains('diff')) %>%
    summary() %>%
    print()

print("Isolation days p_80")
dat.iso.80 %>%
    select(contains('diff')) %>%
    summary() %>%
    print()

# ---- vax prevention -----

# Number of cases prevented by one vaccine:
cpv <- function(df) {
    df$cpv <- df$diff.n.s / df$n.vax
    return(df)
}

vp.40 <- read.csv('out-vax-prev-one-frailty-40-prm_vax_p_0.csv') %>%
    cpv()

vp.80 <- read.csv('out-vax-prev-one-frailty-80-prm_vax_p_0.csv') %>%
    cpv()

cpv.40 <- summary(vp.40$cpv)
cpv.80 <- summary(vp.80$cpv)

print("Mean # of cases averted by 1 vacc dose (p_40):")
print(cpv.40)

print("Mean # of cases averted by 1 vacc dose (p_80):")
print(cpv.80)


print("Mean # of vacc required to avert 1 sympt case (p_40):")
print(1/cpv.40)

print("Mean # of vacc required to avert 1 sympt case (p_80):")
print(1/cpv.80)

