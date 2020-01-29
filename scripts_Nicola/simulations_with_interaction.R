


ph$st.edu_male=ph$eduYears
ph$st.edu_male[which(!is.na(ph$st.edu_new2_male) & ph$dummy.sex==0)]=0
ph$st.edu_female=ph$eduYears
ph$st.edu_female[which(!is.na(ph$st.edu_new2_female) & ph$dummy.sex==1)]=0

ph$st.bmi_male=ph$bmi
ph$st.bmi_male[which(!is.na(ph$st.bmi_male) & ph$dummy.sex==0)]=0
ph$st.bmi_female=ph$bmi
ph$st.bmi_female[which(!is.na(ph$st.bmi_female) & ph$dummy.sex==1)]=0




OR.edu=4 ## Bias OR for edu
OR.bmi=4 ## Bias OR for bmi

set.seed(123456)
z.pmen=(ph$eduYears*log(OR.edu))*ph$dummy.sex+(ph$bmi*log(OR.bmi))*ph$dummy.sex
prob.men=1/(1+exp(-z.pmen))
ph$select=rbinom(n=nrow(ph),size = 1,prob = prob.men)

ph$st.edu_new1=ph$eduYears
ph$st.bmi_new1=ph$bmi
ph$sex_new1=ph$dummy.sex

ph$st.edu_new1[ph$select==0]=NA
ph$st.bmi_new1[ph$select==0]=NA
ph$sex_new1[ph$select==0]=NA
#### Create interaction term

ph$st.edu_new1_male=ph$st.edu_new1
ph$st.edu_new1_male[which(!is.na(ph$st.edu_new1_male) & ph$dummy.sex==0)]=0
ph$st.edu_new1_female=ph$st.edu_new1
ph$st.edu_new1_female[which(!is.na(ph$st.edu_new1_female) & ph$dummy.sex==1)]=0

ph$st.bmi_new1_male=ph$st.bmi_new1
ph$st.bmi_new1_male[which(!is.na(ph$st.bmi_new1_male) & ph$dummy.sex==0)]=0
ph$st.bmi_new1_female=ph$st.bmi_new1
ph$st.bmi_new1_female[which(!is.na(ph$st.bmi_new1_female) & ph$dummy.sex==1)]=0



## Invert man and women

OR.edu=0.25 ## Bias OR for edu
OR.bmi=0.25 ## Bias OR for bmi

set.seed(123456)
z.pmen=(ph$eduYears*log(OR.edu))*ph$dummy.sex+(ph$bmi*log(OR.bmi))*ph$dummy.sex
prob.men=1/(1+exp(-z.pmen))
ph$select=rbinom(n=nrow(ph),size = 1,prob = prob.men)

ph$st.edu_new2=ph$eduYears
ph$st.bmi_new2=ph$bmi
ph$sex_new2=ph$dummy.sex

ph$st.edu_new2[ph$select==0]=NA
ph$st.bmi_new2[ph$select==0]=NA
ph$sex_new2[ph$select==0]=NA
#### Create interaction term

ph$st.edu_new2_male=ph$st.edu_new2
ph$st.edu_new2_male[which(!is.na(ph$st.edu_new2_male) & ph$dummy.sex==0)]=0
ph$st.edu_new2_female=ph$st.edu_new2
ph$st.edu_new2_female[which(!is.na(ph$st.edu_new2_female) & ph$dummy.sex==1)]=0

ph$st.bmi_new2_male=ph$st.bmi_new2
ph$st.bmi_new2_male[which(!is.na(ph$st.bmi_new2_male) & ph$dummy.sex==0)]=0
ph$st.bmi_new2_female=ph$st.bmi_new2
ph$st.bmi_new2_female[which(!is.na(ph$st.bmi_new2_female) & ph$dummy.sex==1)]=0



### Trait list to analyse
# eduYears
# st.bmi_male
# st.bmi_female
# st.edu_male
# st.edu_female
# bmi
# dummy.sex
# st.edu_new1
# st.bmi_new1
# sex_new1
# sex_new2
# st.edu_new1_male
# st.edu_new1_female
# st.bmi_new1_male
# st.bmi_new1_female
#
# st.edu_new2
# st.bmi_new2
# sex_new2
# st.edu_new2_male
# st.edu_new2_female
# st.bmi_new2_male
# st.bmi_new2_female





