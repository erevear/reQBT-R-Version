  library(shiny)
  frequencies<-read.csv("Allele_Frequencies.csv",header = FALSE, stringsAsFactors = FALSE)
  drop_out_frequencies<-read.csv("Drop_Out_Rates.csv", header = FALSE, stringsAsFactors = FALSE)
  #drop_in_frequencies<-read.csv("Drop_In_Rates.csv", header = FALSE, stringsAsFactors = FALSE)
  #drop_in_frequencies<-c(PC0 = .975, PC1 = .02, PC2 = .005)
  colnames(drop_out_frequencies)<-c("Type", "Locus", 25, 50, 100, 150, 250, 500, 6.25, 12.5)
  amounts<-c(NA, NA, 25,50,100,150,250,500, 6.25, 12.5)
  
  
  
  check_dup <- function(x,y){
    #print("x from check dup")
    #print(x) 
    #print("y from check dup")
    #print(y)
    z = (as.numeric(x)^.675)+(as.numeric(y)^.675)
    
  }
  
  
  getFreq<-function(locus, alleles) {
    y<-frequencies[frequencies$V1 == locus & frequencies$V2 %in% alleles, ]
    #y<-y[y$V2 == alleles,] 
    
    return(y)
    
  }
  
  get_allele_freq<-function(locus, race) {
    allele_freq<-locus[,race]
    return(allele_freq)
  }
  
  get_allele_freq_2<-function(alleles, race) {
    freq<-alleles[, race]
    temp<-as.numeric(alleles[,"V2"])
    print("temp")
    print(temp)
    the_names<-c(temp, 3.3)
    print(the_names[1])
    print("the names first occurence")
    print(the_names)
    w<-1-sum(freq)
    freq[length(freq)+1]<-w
    names(freq)<-the_names
    return(freq)
  }
  
  match_mark<-function(combos, alleles) {
    the_names<-names(alleles)
    temp1<-combos[,"allele1"]
    temp2<-combos[,"allele2"]
    #print("names")
    print("alleles from match_mark")
    print(alleles)
    #print(the_names[1])
    for(n in 1:length(alleles)) {
      temp1<-replace(as.numeric(temp1), as.numeric(temp1) == as.numeric(the_names[n]), alleles[n])
      temp2<-replace(as.numeric(temp2), as.numeric(temp2) == as.numeric(the_names[n]), alleles[n])
    }
    combos[,"allele1 freq"]<-temp1
    combos[,"allele2 freq"]<-temp2
    return(combos)
  }
  
  check_het<-function(x, y){
    hom<-0
    if(x==y){
      hom<-1
    }
    return(hom)
  }
  
  generate_allele_combos<-function(alleles, the_names) {
    print("alleles from generate combos to check names")
    print(the_names)
    #names(alleles)<-c(1:length(alleles))
    #print(alleles[1])
    combos<-expand.grid(allele1 = the_names, allele2 = the_names, mark = 0)
    print("combos before mark")
    print(combos)
    combos[,3] <- mapply(check_dup, combos[,1], combos[,2])
    combos[,4]<-mapply(check_het, combos[,1], combos[,2])
    #print("combos after mark")
    #print(combos)
    
    filter<-duplicated(combos[,3])
    
    remove<-which(filter==TRUE)
    #the_names<-names(alleles)
    #print("names")
    
    #print(the_names[1])
    combos<-combos[-c(remove),]
    combos<-match_mark(combos, alleles)
    #for(n in 1:length(alleles)) {
     # combos[,"allele1"]<-replace(as.numeric(combos[,"allele1"]), as.numeric(combos[,"allele1"]) == as.numeric(the_names[n]), alleles[n])
      #combos[,"allele2"]<-replace(as.numeric(combos[,"allele2"]), as.numeric(combos[,"allele2"]) == as.numeric(the_names[n]), alleles[n])
    #}
    #print(new1)
    #print(new2)
    
    return(combos)
  }
  
  generate_contributor_combos<-function(allele_combos, num_contributors) {
    #temp<-paste(allele_combos[,"allele1"], allele_combos[,"allele2"], "")
    if(num_contributors == 1) {
      
      all_combos<-expand.grid(u1=allele_combos[,"mark"])
      #all_combos<-expand.grid(u1=temp)
    }else if(num_contributors == 2) {
      all_combos<-expand.grid(u2=allele_combos[,"mark"], u1=allele_combos[,"mark"])
      #all_combos<-expand.grid(u2=temp, u1=temp)
    }else{
      all_combos<-expand.grid(u3=allele_combos[,"mark"], u2=allele_combos[,"mark"], u1=allele_combos[,"mark"])
      #all_combos<-expand.grid(u3=temp, u2=temp, u1=temp)
    }
    
    return(all_combos)
  }
  
  
  get_dropout_range<-function(num_contributors, D_ND) {
    Het1_loc<-paste("HET1", num_contributors, D_ND, sep="-")
    Het2_loc<-paste("HET2", num_contributors, D_ND, sep="-")
    Hom1_loc<-paste("HOM1", num_contributors, D_ND, sep="-")
    Het1_table<-drop_out_frequencies[drop_out_frequencies$Type == Het1_loc,]
    #Het1_table[16,]<-amounts
    Het2_table<-drop_out_frequencies[drop_out_frequencies$Type == Het2_loc,]
    #Het2_table[16,]<-amounts
    Hom1_table<-drop_out_frequencies[drop_out_frequencies$Type == Hom1_loc,]
    #Hom1_table[16,]<-amounts
    #colnames(Hom1_table)<-c("Type", "Locus.1", 25, 50, 100, 150, 250, 500, 6.25, 12.5)
    #colnames(Het1_table)<-c("Type", "Locus.2", 25, 50, 100, 150, 250, 500, 6.25, 12.5)
    #colnames(Het2_table)<-c("Type", "Locus.3", 25, 50, 100, 150, 250, 500, 6.25, 12.5)
    all_temp<-cbind(Hom1_table, Het1_table, Het2_table)
    
    #print("all temp")
    #print(all_temp)
    return(all_temp)
    #drop_out_frequencies[drop_out_frequencies$Type == "HET1-4-ND",]
  }
  
  calculate_drop_out<-function(dropout_rate_table, quant, locus) {
    #if(quant > 500){
     # quant <- 500
    #}
    #temp_table1<-dropout_rate_table[dropout_rate_table$Locus.1 == locus,]
    #temp_table2<-dropout_rate_table[dropout_rate_table$Locus.2 == locus,]
    #temp_table3<-dropout_rate_table[dropout_rate_table$Locus.3 == locus,]
    temp_table<-dropout_rate_table[dropout_rate_table$Locus == locus,]
    print("temp_table")
    print(temp_table)
    print("quant from dropout calculate")
    print(quant)
    possible_amounts<-c( 6.25, 12.5, 25,50,100,150,250,500)
    hom1<-temp_table[,1:10]
    het1<-temp_table[,11:20]
    het2<-temp_table[,21:30]
    #temp_table[2,]<-amounts
    if(quant%in%possible_amounts){
      Hom1_do_freq<-as.numeric(hom1[,toString(quant)])
      Het1_do_freq<-as.numeric(het1[,toString(quant)])
      Het2_do_freq<-as.numeric(het2[,toString(quant)])
    }else{
      
      possible_amounts[length(possible_amounts)+1]<-as.numeric(quant)
      possible_amounts<-sort(possible_amounts)
      quant_left<-possible_amounts[which(possible_amounts== as.numeric(quant)) - 1]
      quant_right<-possible_amounts[which(possible_amounts == as.numeric(quant)) + 1]
      hom1_quant_right<-as.numeric(hom1[,toString(quant_right)])
      hom1_quant_left<-as.numeric(hom1[,toString(quant_left)])
    
      het1_quant_right<-as.numeric(het1[,toString(quant_right)])
      het1_quant_left<-as.numeric(het1[,toString(quant_left)])
    
      het2_quant_right<-as.numeric(het2[,toString(quant_right)])
      het2_quant_left<-as.numeric(het2[,toString(quant_left)])
      quant<-as.numeric(quant)
      hom1_slope<-(hom1_quant_right - hom1_quant_left)/(quant_right-quant_left)
      het1_slope<-(het1_quant_right - het1_quant_left)/(quant_right-quant_left)
      het2_slope<-(het2_quant_right - het2_quant_left)/(quant_right-quant_left)
      
      hom1_b_intercept<-hom1_quant_left - hom1_slope * quant_left
      het1_b_intercept<-het1_quant_left - het1_slope * quant_left
      het2_b_intercept<-het2_quant_left - het2_slope * quant_left
      
      Hom1_do_freq<-hom1_slope * quant + hom1_b_intercept
      Het1_do_freq<-het1_slope * quant + het1_b_intercept
      Het2_do_freq<-het2_slope * quant + het2_b_intercept
      #Hom1_do_freq<-hom1_quant_right+(hom1_quant_left - hom1_quant_right) * (as.numeric(quant_right) - as.numeric(quant))/(as.numeric(quant_left)-as.numeric(quant_right))
      #Het1_do_freq<-het1_quant_right+(het1_quant_left - het1_quant_right) * (as.numeric(quant_right) - as.numeric(quant))/(as.numeric(quant_left)-as.numeric(quant_right))
      #Het2_do_freq<-het2_quant_right+(het2_quant_left - het2_quant_right) * (as.numeric(quant_right) - as.numeric(quant))/(as.numeric(quant_left)-as.numeric(quant_right))
    #Hom1 = ND_3_Hom1[allele_match,quant_right]+(ND_3_Hom1[allele_match,quant_left]-ND_3_Hom1[allele_match,quant_right])*(ND_3_Hom1[16,quant_right]-quant)/(ND_3_Hom1[16,quant_left]-ND_3_Hom1[16,quant_right])
    }
    Hom0_do_freq<- 1-Hom1_do_freq
    Het0_do_freq<- 1-(Het1_do_freq + Het2_do_freq)
    
    Hom_dropouts<-c(Hom1_do_freq,0, Hom0_do_freq)
    Het_dropouts<-c(Het2_do_freq, Het1_do_freq, Het0_do_freq)
    dropouts<-rbind(Hom_dropouts, Het_dropouts)
    
    #print(hom1_quant_left)
    #print(hom1_quant_right)
    #print(hom1)
    return(dropouts)
    
    
    

  }
  
  generate_freq_table<-function(allele_combos) {
    allele_combos[,"probability"] <- mapply(geno_freq, as.numeric(allele_combos[,"allele1"]), as.numeric(allele_combos[,"allele2"]), as.numeric(allele_combos[,"allele1 freq"]), as.numeric(allele_combos[,"allele2 freq"]))
    return(allele_combos)
  }
  
  geno_freq <- function(x,y,a,b){
    
    if(x==y){
      
      z<- a^2+0.03*a*(1-a)
    }
    
    else{
      z<- 2*a*b
    }
  }
  
  generate_rep_dropout_table<-function(allele_combos, dropout_rates, rep_info) {
    rates<-rep(NA, ncol(allele_combos))
    for(n in 1:nrow(allele_combos)) {
      if(as.numeric(allele_combos[n,"allele1"] == allele_combos[n,"allele2"])){
      #if(as.numeric(allele_combos[n, 1]) == as.numeric(allele_combos[n, 2])) {
        allele_combos[n,"dropout"]<-dropout_rates[1, (1+(allele_combos[n,"allele1"]%in%rep_info)+allele_combos[n,"allele2"]%in%rep_info)]
      }else{
        allele_combos[n,"dropout"]<-dropout_rates[2, (1+(allele_combos[n,"allele1"]%in%rep_info)+allele_combos[n,"allele2"]%in%rep_info)]
      }
    }
    return(allele_combos)
  }
  
  get_full_freq_table<-function(allele_combos, all_contributor_combos) {
    for(n in 1:nrow(allele_combos)){
      for(m in 1:ncol(all_contributor_combos)){
        all_contributor_combos[,m]<-replace(all_contributor_combos[,m], all_contributor_combos[,m]==allele_combos[n,"mark"], allele_combos[n,"probability"])
      }
    }
    return(all_contributor_combos)
  }
  
  get_full_rep_table<-function(allele_combos, all_contributor_combos){
    for(n in 1:nrow(allele_combos)){
      for(m in 1:ncol(all_contributor_combos)){
        all_contributor_combos[,m]<-replace(all_contributor_combos[,m], all_contributor_combos[,m]==allele_combos[n,"mark"], allele_combos[n,"dropout"])
      }
    }
    return(all_contributor_combos)
  }
  
  get_full_allele_table<-function(allele_combos, all_contributor_combos){
      for(m in 1:ncol(all_contributor_combos)){
        for(n in 1:nrow(all_contributor_combos)){
          all_contributor_combos[n, paste(m,"allele1")]<-allele_combos[allele_combos$mark==all_contributor_combos[n,m], 1]
          all_contributor_combos[n, paste(m,"allele2")]<-allele_combos[allele_combos$mark==all_contributor_combos[n,m], 2]
        }
      }
      return(all_contributor_combos)
  }
  
  calculate_numerator<-function(knowns, full_table, num_cont, allele_combos, names){
   
   conts<-rep(NA,(length(knowns)/2))
   print("conts empty")
   print(conts)
   l<-1
   print("length knowns")
   print(length(knowns))
   for(n in 1:length(conts)){
     conts[n]<-check_dup(as.numeric(names[l]), as.numeric(names[l+1]))
       #get_mark_with_alleles(allele_combos, knowns[l:(l+1)])
     l<-l+2
     #n<-n+1
   }
   print("contributors")
   print(conts)
   #print(known_contributor1[,"V3"])
   #print("contributor 2")
   #print(known_contributor2[,"V3"])
   #cont1<-get_mark_with_alleles(allele_combos, known_contributor1[,"V3"])
   #cont2<-get_mark_with_alleles(allele_combos, known_contributor2[,"V3"])
   #print("cont1")
   #print(cont1)
   #print("cont2")
   #print(cont2)
   #conts<-c(cont1, cont2)
   m<-num_cont
   for(n in 1:length(conts)){
      full_table<-full_table[full_table[,m]==conts[n],]
      m<-m-1
   }
   print("table before from numerator multiplication")
   print(full_table[1:15,])
   full_table<-full_table[,(num_cont + (num_cont*2)+1):ncol(full_table)]
   temp1<-full_table[,1]
   temp2<-full_table[,num_cont]
   full_table[,1]<-temp2
   full_table[,num_cont]<-temp1
   num<-length(conts)+1
   print("num")
   print(num)
   col<-ncol(full_table)+1
   print("col")
   print(col)
   for(n in 1:nrow(full_table)){
     full_table[n,col]<-prod(full_table[n, num:(col-1)])
   }
   print("full table after multiplication")
   print(full_table[1:15,])
   print("probability")
   return(sum(full_table[,ncol(full_table)]))
    
  }
  
  calculate_denominator<-function(knowns, allele_combos, full_table, num_cont, names){
    conts<-0
    num<-0
    if(length(knowns)>0){
      conts<-rep(NA,(length(knowns)/2))
      print("conts empty")
      print(conts)
      l<-1
      print("length knowns")
      print(length(knowns))
      for(n in 1:length(conts)){
        conts[n]<-check_dup(as.numeric(names[l]), as.numeric(names[l+1]))
          #get_mark_with_alleles(allele_combos, knowns[l:(l+1)])
        l<-l+2
      #n<-n+1
      }
      print("contributors")
      print(conts)
    #print(known_contributor1[,"V3"])
    #print("contributor 2")
    #print(known_contributor2[,"V3"])
    #cont1<-get_mark_with_alleles(allele_combos, known_contributor1[,"V3"])
    #cont2<-get_mark_with_alleles(allele_combos, known_contributor2[,"V3"])
    #print("cont1")
    #print(cont1)
    #print("cont2")
    #print(cont2)
    #conts<-c(cont1, cont2)
      m<-num_cont
      for(n in 1:length(conts)){
        full_table<-full_table[full_table[,m]==conts[n],]
        m<-m-1
      }
    num<-length(conts)+1
    }
    print("table before multiplication from denominator")
    full_table<-full_table[,(num_cont + (num_cont*2)+1):ncol(full_table)]
    temp1<-full_table[,1]
    temp2<-full_table[,num_cont]
    full_table[,1]<-temp2
    full_table[,num_cont]<-temp1
    
    col<-ncol(full_table)+1
    for(n in 1:nrow(full_table)){
      full_table[n,col]<-prod(full_table[n, num:(col-1)])
    }
    print("full table from denominator calculation")
    print(full_table[1:15,])
    print("probability denominator")
    return(sum(full_table[,ncol(full_table)]))
  }
  
  get_mark_with_alleles<-function(allele_combos, alleles) {
    print("alleles from get mark")
    print(alleles)
    row<-allele_combos[(allele_combos$allele1 == alleles[1] & allele_combos$allele2 == alleles[2]) | (allele_combos$allele2 == alleles[1] & allele_combos$allele1 == alleles[2]) ,]
    print("allele combos from get mark")
    print(allele_combos)
    print("row")
    print(row)
    
    return(row["mark"])
  }
  
  generate_full_dropout<-function(rep_1, rep_2, rep_3, allele_combos, dropout_freq_table, all_the_combos){
   full_table<-rep(NA, nrow(all_the_combos))
    if(!("INC"%in%rep_1) & length(rep_1) > 0){
      Rep_1_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, rep_1)
      print("dropout tables")
      #print(Rep_1_dropout_table)
      full_rep1_table<-get_full_rep_table(Rep_1_dropout_table, all_the_combos)
      print(full_rep1_table)
      full_table<-cbind(full_table, full_rep1_table)
      #full_table<-full_rep1_table
    }
    if(!("INC"%in%rep_2) & length(rep_2) > 0){
      Rep_2_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, rep_2)
      full_rep2_table<-get_full_rep_table(Rep_2_dropout_table, all_the_combos)
      print(full_rep2_table)
      full_table<-cbind(full_table, full_rep2_table)
    }
    if(!("INC"%in%rep_3) & length(rep_3) > 0){
      Rep_3_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, rep_3)
      full_rep3_table<-get_full_rep_table(Rep_3_dropout_table, all_the_combos)
      print(full_rep3_table)
      full_table<-cbind(full_table, full_rep3_table)
    }
    #full_table<-cbind(full_rep1_table, full_rep2_table, full_rep3_table)
    full_table<-subset(full_table, select=-full_table)
    print(full_table)
    return(full_table)
  }
  
  generate_full_dropin_table<-function(rep_1, rep_2, rep_3, full_allele_table, number_contributors, deducible_nondeducible){
    full_table<-rep(NA, nrow(full_allele_table))
    if(!("INC"%in%rep_1) & length(rep_1) > 0){
      
      full_rep1_table<-apply_dropin(full_allele_table, rep_1, number_contributors, deducible_nondeducible)
      full_table<-cbind(full_table, full_rep1_table[,ncol(full_rep1_table)])
      #full_table<-full_rep1_table
    }
    if(!("INC"%in%rep_2) & length(rep_2) > 0){
      full_rep2_table<-apply_dropin(full_allele_table, rep_2, number_contributors, deducible_nondeducible)
      full_table<-cbind(full_table, full_rep2_table[,ncol(full_rep2_table)])
    }
    if(!("INC"%in%rep_3) & length(rep_3) > 0){
      full_rep3_table<-apply_dropin(full_allele_table, rep_3, number_contributors, deducible_nondeducible)
      full_table<-cbind(full_table, full_rep3_table[,ncol(full_rep3_table)])
    }
    #full_table<-cbind(full_rep1_table, full_rep2_table, full_rep3_table)
    full_table<-subset(full_table, select=-full_table)
    #print(full_table)
    return(full_table)
  }
  
  
  make_full_locus_table<-function(locus_alleles,rep_1_alleles, rep_2_alleles, rep_3_alleles, locus, number_contributors, deducible_nondeducible, quantity, n_known1, n_known2, n_known3, 
                                  d_known1, d_known2, d_known3, race){
    dropout_table<-get_dropout_range(number_contributors, deducible_nondeducible)
    #print("rep 1 alleles")
    #print(rep_1_alleles)
    print("dropout table")
    print(dropout_table)
    dropout_freq_table<-calculate_drop_out(dropout_table, quantity, locus)
    #print("locus dropout calculated")
    print(dropout_freq_table)
    #allele_frequency_rep1<-get_allele_freq_2(rep_1_alleles, race)
    #allele_frequency_rep2<-get_allele_freq_2(rep_2_alleles, race)
    #allele_frequency_rep3<-get_allele_freq_2(rep_3_alleles, race)
    #print("the allele frequency")
   # print(allele_frequency_rep1)
    #print(allele_frequency_rep2)
    #print(allele_frequency_rep3)
    present_allele_freq<-get_allele_freq_2(locus_alleles, race)
    print("present alleles")
    print(present_allele_freq)
    the_names<-names(present_allele_freq)
    allele_combos<-generate_allele_combos(present_allele_freq, as.numeric(the_names))
    #allele_combos_rep1<-generate_allele_combos(allele_frequency_rep1)
    #allele_combos_rep2<-generate_allele_combos(allele_frequency_rep2)
    #allele_combos_rep3<-generate_allele_combos(allele_frequency_rep3)
    print("allele combos")
    print(allele_combos)
  pg_table<-generate_freq_table(allele_combos)
    #pg_table_rep1<-generate_freq_table(allele_combos_rep1)
    #pg_table_rep2<-generate_freq_table(allele_combos_rep2)
    #pg_table_rep3<-generate_freq_table(allele_combos_rep3)
  print("pg table")
  print(pg_table)
  all_the_combos<-generate_contributor_combos(allele_combos, number_contributors)
  #print("contributor combos")
  #print(all_the_combos)
    #print("full freq table")
  full_frequency_table<-get_full_freq_table(pg_table, all_the_combos)
  print(full_frequency_table[1:100,])
  print("rep 1 allles")
  print(rep_1_alleles)
  #Rep_1_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, rep_1_alleles[,"V2"])
  #print("rep one drop table")
  #print(Rep_1_dropout_table)
  #Rep_2_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, rep_2_alleles[,"V2"])
  #Rep_3_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, rep_3_alleles[,"V2"])
    #print("rep 1 dropout table")
    #print(Rep_1_dropout_table)
  #full_rep1_table<-get_full_rep_table(Rep_1_dropout_table, all_the_combos)
    #print("rep1 do table")
    #print(full_rep1_table)
  #full_rep2_table<-get_full_rep_table(Rep_2_dropout_table, all_the_combos)
  #full_rep3_table<-get_full_rep_table(Rep_3_dropout_table, all_the_combos)
    #print("full rep table")
    #print(full_rep_table)
  full_dropout_table<-generate_full_dropout(rep_1_alleles[,"V2"], rep_2_alleles[,"V2"], rep_3_alleles[,"V2"], allele_combos, dropout_freq_table, all_the_combos)
  print("full dropout table")
  print(full_dropout_table)
  full_allele_table<-get_full_allele_table(allele_combos, all_the_combos)
  full_dropin_table<-generate_full_dropin_table(rep_1_alleles[,"V2"], rep_2_alleles[,"V2"], rep_3_alleles[,"V2"], full_allele_table, number_contributors, deducible_nondeducible)
  print("full dropin table")
  print(full_dropin_table)
 
  #print("full allele table")
  #print(full_allele_table[1:100,])
  #rep1_with_drop_in<-apply_dropin(full_allele_table, rep_1_alleles[,"V2"], number_contributors, deducible_nondeducible)
  #print("rep one with drop in")
  #print(rep1_with_drop_in)
  #rep2_with_drop_in<-apply_dropin(full_allele_table, rep_2_alleles[,"V2"], number_contributors, deducible_nondeducible)
  #rep3_with_drop_in<-apply_dropin(full_allele_table, rep_3_alleles[,"V2"], number_contributors, deducible_nondeducible)
    #print("drop in")
    #print(with_drop_in)
  everything<-cbind(full_allele_table, full_frequency_table, full_dropin_table, full_dropout_table)
  print(everything[1:150,])
    #print("knowns")
    #print(n_known1)
    #print(n_known2)
    #print(n_known3)
  n_knowns<-c(check_knowns(n_known1[,race], allele_combos), check_knowns(n_known2[,race],allele_combos), check_knowns(n_known3[,race],allele_combos))
  n_knowns_names<-c(check_knowns(n_known1[,"V2"],allele_combos), check_knowns(n_known2[,"V2"],allele_combos), check_knowns(n_known3[,"V2"],allele_combos))
   
    #print("numerator knowns V6")
    #print(n_knowns)
    #knowns<-knowns[!is.na(knowns)]
    
    #print("denominator knowns")
    #print(d_known1)
    #print(d_known2)
    #print(d_known3)
  d_knowns<-c(check_knowns(d_known1[,race],allele_combos), check_knowns(d_known2[,race],allele_combos), check_knowns(d_known3[,race], allele_combos))
  d_knowns_names<-c(check_knowns(d_known1[,"V2"],allele_combos), check_knowns(d_known2[,"V2"],allele_combos), check_knowns(d_known3[,"V2"],allele_combos))
    #print("denominator knowns V6")
    #print(d_knowns)
    
  numerator<-calculate_numerator(n_knowns, everything, number_contributors, allele_combos, n_knowns_names)
  denominator<-calculate_denominator(d_knowns, allele_combos, everything, number_contributors, d_knowns_names)
  
    #print(calculate_numerator(known1, known2, everything, number_contributors, allele_combos))
    #print(calculate_denominator(everything, number_contributors))
  print("numerator prob")
  print(numerator)
  print("denominator prob")
  print(denominator)
  LR<-numerator/denominator
  print("LR")
  print(LR)
  #return(LR)
  }
  
  check_knowns<-function(known, allele_combos){
    if(length(known) == 1){
      known[2]<-known[1]
    }
    if(length(known) >= 1){
      for(n in 1: length(known)){
        if(!(known[n]%in%allele_combos[,"allele1"]) & !(known[n]%in%allele_combos[,"allele2"])){
          known[n]<-3.3
        }
      }
    }
    return(known)
  }
  
  generate_name<-function(known_names){
    if(length(known_names) == 1){
      known_names[2]<-known_names[1]
    }
    return(known_names)
  }
  
  apply_dropin<-function(full_allele_table, rep, num_contributors, D_ND) {
    
    if(D_ND == "ND"){
      drop_in_frequencies<-c(PC0 = .975, PC1 = .02, PC2 = .005)
    }else{
      drop_in_frequencies<-c(PC0 = .96, PC1 = .035, PC2 = .005)
    }
    columns<-ncol(full_allele_table)
    print("rep from drop in part")
    print(rep)
    print("drop in from apply_dropin")
    print(drop_in_frequencies)
    for(n in 1:nrow(full_allele_table)){
      alleles<-as.numeric(full_allele_table[n,(num_contributors + 1) :columns])
      #print("alleles from drop in part")
      #print(alleles)
      rep_in_alleles<-(rep%in%alleles)
      
      count<-length(which(rep_in_alleles==TRUE))
      i<-1 + (length(rep) - count)
      #print("rep in alleles")
      #print(rep%in%alleles)
      #print("i")
      #print(i)
      if(i <= 3) {
        drop_in<-drop_in_frequencies[i]
      }else{
        drop_in<-drop_in_frequencies[3]
      }
      full_allele_table[n, columns +1]<-drop_in
      
      
    }
    return(full_allele_table)
  }
  

  
  ui<-shinyUI(fluidPage(
                
                mainPanel(
                  tabsetPanel(
                    tabPanel("Case Information", textInput(inputId = "Quantity", label = "Enter the quantity of DNA"), 
                                 numericInput(inputId = "numCont", label = "Select the number of contributors (1 - 4)", value = 1, min = 1, max = 4, step = NA),
                                 selectInput(inputId = "d_nd", label = "Is this profile deducible or non-deducible", choices = c("Deducible" = "D", "Non-Deducible" = "ND"), selected = "D", multiple = FALSE),
                                 actionButton(inputId = "submitButton", label = "Submit")),
                tabPanel("Denominator", tags$h3("Enter alleles of known contributors in the denominator"), tags$h4("U1"), checkboxGroupInput(inputId = "d_known_u_one_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D7", label = "D7",choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_vWA", label = "vWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                
                checkboxGroupInput(inputId = "d_known_u_one_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_one_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
                tags$h4("U2"),
                checkboxGroupInput(inputId = "d_known_u_two_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_vWA", label = "vWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                
                checkboxGroupInput(inputId = "d_known_u_two_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_two_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
                tags$h4("U3"),
                checkboxGroupInput(inputId = "d_known_u_three_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_vWA", label = "vWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                
                checkboxGroupInput(inputId = "d_known_u_three_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                checkboxGroupInput(inputId = "d_known_u_three_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE)
                #actionButton(inputId = "submitButton", label = "Submit")
                
                ),
                
                tabPanel("Numerator", tags$h3("Enter alleles of known contributors in the numerator"), tags$h4("U1"), checkboxGroupInput(inputId = "n_known_u_one_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                         
                         checkboxGroupInput(inputId = "n_known_u_one_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_one_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
                         tags$h4("U2"),
                         checkboxGroupInput(inputId = "n_known_u_two_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                         
                         checkboxGroupInput(inputId = "n_known_u_two_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_two_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
                         tags$h4("U3"),
                         checkboxGroupInput(inputId = "n_known_u_three_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_three_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_three_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_three_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_FGA", label = "FGA",choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                         
                         checkboxGroupInput(inputId = "n_known_u_three_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "n_known_u_three_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE)
                         #actionButton(inputId = "submitButton", label = "Submit")
                         
                ),
                
                tabPanel("Replicates", tags$h3("Enter replicate results"),tags$h4("Rep 1"), checkboxGroupInput(inputId = "rep_one_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D21", label = "D21",choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                         
                         checkboxGroupInput(inputId = "rep_one_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_one_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
                         tags$h4("Rep 2"),
                         checkboxGroupInput(inputId = "rep_two_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D8", label = "D8",choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                         
                         checkboxGroupInput(inputId = "rep_two_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_two_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
                         tags$h4("Rep 3"),
                         checkboxGroupInput(inputId = "rep_three_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
                         
                         checkboxGroupInput(inputId = "rep_three_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
                         checkboxGroupInput(inputId = "rep_three_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE)
                         
                         
                )
                
                ))
                
                #actionButton(inputId = "submitButton", label = "Submit")
                ))
  
  
  server<- function(input, output){
    observeEvent(input$submitButton, {
      quantity<-input$Quantity
      number_contributors<-input$numCont
      deducible_nondeducible<-input$d_nd
      
      #print("location in do table")
      #dropout_table<-get_dropout_range(number_contributors, deducible_nondeducible)
      #print(dropout_table)
      #dropout_freq_table<-calculate_drop_out(dropout_table, quantity, "TPOX")
      #print(dropout_freq_table)
    
      #print(quantity)
      #print(number_contributors)
      #print(deducible_nondeducible)
      d_TPOX_u_one<-as.numeric(input$d_known_u_one_TPOX)
      d_CSF_u_one<-as.numeric(input$d_known_u_one_CSF)
      d_D3_u_one<-as.numeric(input$d_known_u_one_D3)
      d_D16_u_one<-as.numeric(input$d_known_u_one_D16)
      d_D7_u_one<-as.numeric(input$d_known_u_one_D7)
      d_D8_u_one<-as.numeric(input$d_known_u_one_D8)
      d_D21_u_one<-as.numeric(input$d_known_u_one_D21)
      d_D18_u_one<-as.numeric(input$d_known_u_one_D18)
      d_FGA_u_one<-as.numeric(input$d_known_u_one_FGA)
      d_D5_u_one<-as.numeric(input$d_known_u_one_D5)
      d_D13_u_one<-as.numeric(input$d_known_u_one_D13)
      d_vWA_u_one<-as.numeric(input$d_known_u_one_vWA)
      d_TH01_u_one<-as.numeric(input$d_known_u_one_TH01)
      d_D2_u_one<-as.numeric(input$d_known_u_one_D2)
      d_D19_u_one<-as.numeric(input$d_known_u_one_D19)
      
      
      d_TPOX_u_two<-as.numeric(input$d_known_u_two_TPOX)
      d_CSF_u_two<-as.numeric(input$d_known_u_two_CSF)
      d_D3_u_two<-as.numeric(input$d_known_u_two_D3)
      d_D16_u_two<-as.numeric(input$d_known_u_two_D16)
      d_D7_u_two<-as.numeric(input$d_known_u_two_D7)
      d_D8_u_two<-as.numeric(input$d_known_u_two_D8)
      d_D21_u_two<-as.numeric(input$d_known_u_two_D21)
      d_D18_u_two<-as.numeric(input$d_known_u_two_D18)
      d_FGA_u_two<-as.numeric(input$d_known_u_two_FGA)
      d_D5_u_two<-as.numeric(input$d_known_u_two_D5)
      d_D13_u_two<-as.numeric(input$d_known_u_two_D13)
      d_vWA_u_two<-as.numeric(input$d_known_u_two_vWA)
      d_TH01_u_two<-as.numeric(input$d_known_u_two_TH01)
      d_D2_u_two<-as.numeric(input$d_known_u_two_D2)
      d_D19_u_two<-as.numeric(input$d_known_u_two_D19)
      
      
      d_TPOX_u_three<-as.numeric(input$d_known_u_three_TPOX)
      d_CSF_u_three<-as.numeric(input$d_known_u_three_CSF)
      d_D3_u_three<-as.numeric(input$d_known_u_three_D3)
      d_D16_u_three<-as.numeric(input$d_known_u_three_D16)
      d_D7_u_three<-as.numeric(input$d_known_u_three_D7)
      d_D8_u_three<-as.numeric(input$d_known_u_three_D8)
      d_D21_u_three<-as.numeric(input$d_known_u_three_D21)
      d_D18_u_three<-as.numeric(input$d_known_u_three_D18)
      d_FGA_u_three<-as.numeric(input$d_known_u_three_FGA)
      d_D5_u_three<-as.numeric(input$d_known_u_three_D5)
      d_D13_u_three<-as.numeric(input$d_known_u_three_D13)
      d_vWA_u_three<-as.numeric(input$d_known_u_three_vWA)
      d_TH01_u_three<-as.numeric(input$d_known_u_three_TH01)
      d_D2_u_three<-as.numeric(input$d_known_u_three_D2)
      d_D19_u_three<-as.numeric(input$d_known_u_three_D19)
      
      
      
      d_TPOX1<-getFreq("TPOX",d_TPOX_u_one)
      d_CSF1<-getFreq("CSF",d_CSF_u_one)
      d_D161<-getFreq("D16",d_D16_u_one)
      d_D51<-getFreq("D5", d_D5_u_one)
      d_vWA1<-getFreq("vWA",d_vWA_u_one)
      d_D31<-getFreq("D3",d_D3_u_one)
      d_D71<-getFreq("D7",d_D7_u_one)
      d_D81<-getFreq("D8",d_D8_u_one)
      d_D181<-getFreq("D18",d_D18_u_one)
      d_FGA1<-getFreq("FGA",d_FGA_u_one)
      d_TH011<-getFreq("TH01",d_TH01_u_one)
      d_D21<-getFreq("D2",d_D2_u_one)
      d_D191<-getFreq("D19",d_D19_u_one)
      d_D211<-getFreq("D21",d_D21_u_one)
      d_D131<-getFreq("D13",d_D13_u_one)
      
      d_TPOX2<-getFreq("TPOX",d_TPOX_u_two)
      d_CSF2<-getFreq("CSF",d_CSF_u_two)
      d_D162<-getFreq("D16",d_D16_u_two)
      d_D52<-getFreq("D5", d_D5_u_two)
      d_vWA2<-getFreq("vWA",d_vWA_u_two)
      d_D32<-getFreq("D3",d_D3_u_two)
      d_D72<-getFreq("D7",d_D7_u_two)
      d_D82<-getFreq("D8",d_D8_u_two)
      d_D182<-getFreq("D18",d_D18_u_two)
      d_FGA2<-getFreq("FGA",d_FGA_u_two)
      d_TH012<-getFreq("TH01",d_TH01_u_two)
      d_D22<-getFreq("D2",d_D2_u_two)
      d_D192<-getFreq("D19",d_D19_u_two)
      d_D212<-getFreq("D21",d_D21_u_two)
      d_D132<-getFreq("D13",d_D13_u_two)
      
      d_TPOX3<-getFreq("TPOX",d_TPOX_u_three)
      d_CSF3<-getFreq("CSF",d_CSF_u_three)
      d_D163<-getFreq("D16",d_D16_u_three)
      d_D53<-getFreq("D5", d_D5_u_three)
      d_vWA3<-getFreq("vWA",d_vWA_u_three)
      d_D33<-getFreq("D3",d_D3_u_three)
      d_D73<-getFreq("D7",d_D7_u_three)
      d_D83<-getFreq("D8",d_D8_u_three)
      d_D183<-getFreq("D18",d_D18_u_three)
      d_FGA3<-getFreq("FGA",d_FGA_u_three)
      d_TH013<-getFreq("TH01",d_TH01_u_three)
      d_D23<-getFreq("D2",d_D2_u_three)
      d_D193<-getFreq("D19",d_D19_u_three)
      d_D213<-getFreq("D21",d_D21_u_three)
      d_D133<-getFreq("D13",d_D13_u_three)
      
      
      n_TPOX_u_one<-as.numeric(input$n_known_u_one_TPOX)
      n_CSF_u_one<-as.numeric(input$n_known_u_one_CSF)
      n_D3_u_one<-as.numeric(input$n_known_u_one_D3)
      n_D16_u_one<-as.numeric(input$n_known_u_one_D16)
      n_D7_u_one<-as.numeric(input$n_known_u_one_D7)
      n_D8_u_one<-as.numeric(input$n_known_u_one_D8)
      n_D21_u_one<-as.numeric(input$n_known_u_one_D21)
      n_D18_u_one<-as.numeric(input$n_known_u_one_D18)
      n_FGA_u_one<-as.numeric(input$n_known_u_one_FGA)
      n_D5_u_one<-as.numeric(input$n_known_u_one_D5)
      n_D13_u_one<-as.numeric(input$n_known_u_one_D13)
      n_vWA_u_one<-as.numeric(input$n_known_u_one_vWA)
      n_TH01_u_one<-as.numeric(input$n_known_u_one_TH01)
      n_D2_u_one<-as.numeric(input$n_known_u_one_D2)
      n_D19_u_one<-as.numeric(input$n_known_u_one_D19)
      
      
      n_TPOX_u_two<-as.numeric(input$n_known_u_two_TPOX)
      n_CSF_u_two<-as.numeric(input$n_known_u_two_CSF)
      n_D3_u_two<-as.numeric(input$n_known_u_two_D3)
      n_D16_u_two<-as.numeric(input$n_known_u_two_D16)
      n_D7_u_two<-as.numeric(input$n_known_u_two_D7)
      n_D8_u_two<-as.numeric(input$n_known_u_two_D8)
      n_D21_u_two<-as.numeric(input$n_known_u_two_D21)
      n_D18_u_two<-as.numeric(input$n_known_u_two_D18)
      n_FGA_u_two<-as.numeric(input$n_known_u_two_FGA)
      n_D5_u_two<-as.numeric(input$n_known_u_two_D5)
      n_D13_u_two<-as.numeric(input$n_known_u_two_D13)
      n_vWA_u_two<-as.numeric(input$n_known_u_two_vWA)
      n_TH01_u_two<-as.numeric(input$n_known_u_two_TH01)
      n_D2_u_two<-as.numeric(input$n_known_u_two_D2)
      n_D19_u_two<-as.numeric(input$n_known_u_two_D19)
      
      
      n_TPOX_u_three<-as.numeric(input$n_known_u_three_TPOX)
      n_CSF_u_three<-as.numeric(input$n_known_u_three_CSF)
      n_D3_u_three<-as.numeric(input$n_known_u_three_D3)
      n_D16_u_three<-as.numeric(input$n_known_u_three_D16)
      n_D7_u_three<-as.numeric(input$n_known_u_three_D7)
      n_D8_u_three<-as.numeric(input$n_known_u_three_D8)
      n_D21_u_three<-as.numeric(input$n_known_u_three_D21)
      n_D18_u_three<-as.numeric(input$n_known_u_three_D18)
      n_FGA_u_three<-as.numeric(input$n_known_u_three_FGA)
      n_D5_u_three<-as.numeric(input$n_known_u_three_D5)
      n_D13_u_three<-as.numeric(input$n_known_u_three_D13)
      n_vWA_u_three<-as.numeric(input$n_known_u_three_vWA)
      n_TH01_u_three<-as.numeric(input$n_known_u_three_TH01)
      n_D2_u_three<-as.numeric(input$n_known_u_three_D2)
      n_D19_u_three<-as.numeric(input$n_known_u_three_D19)
      
      
      
      n_TPOX1<-getFreq("TPOX",n_TPOX_u_one)
      n_CSF1<-getFreq("CSF",n_CSF_u_one)
      n_D161<-getFreq("D16",n_D16_u_one)
      n_D51<-getFreq("D5", n_D5_u_one)
      n_vWA1<-getFreq("vWA",n_vWA_u_one)
      n_D31<-getFreq("D3",n_D3_u_one)
      n_D71<-getFreq("D7",n_D7_u_one)
      n_D81<-getFreq("D8",n_D8_u_one)
      n_D181<-getFreq("D18",n_D18_u_one)
      n_FGA1<-getFreq("FGA",n_FGA_u_one)
      n_TH011<-getFreq("TH01",n_TH01_u_one)
      n_D21<-getFreq("D2",n_D2_u_one)
      n_D191<-getFreq("D19",n_D19_u_one)
      n_D211<-getFreq("D21",n_D21_u_one)
      n_D131<-getFreq("D13",n_D13_u_one)
      
      n_TPOX2<-getFreq("TPOX",n_TPOX_u_two)
      n_CSF2<-getFreq("CSF",n_CSF_u_two)
      n_D162<-getFreq("D16",n_D16_u_two)
      n_D52<-getFreq("D5", n_D5_u_two)
      n_vWA2<-getFreq("vWA",n_vWA_u_two)
      n_D32<-getFreq("D3",n_D3_u_two)
      n_D72<-getFreq("D7",n_D7_u_two)
      n_D82<-getFreq("D8",n_D8_u_two)
      n_D182<-getFreq("D18",n_D18_u_two)
      n_FGA2<-getFreq("FGA",n_FGA_u_two)
      n_TH012<-getFreq("TH01",n_TH01_u_two)
      n_D22<-getFreq("D2",n_D2_u_two)
      n_D192<-getFreq("D19",n_D19_u_two)
      n_D212<-getFreq("D21",n_D21_u_two)
      n_D132<-getFreq("D13",n_D13_u_two)
      
      n_TPOX3<-getFreq("TPOX",n_TPOX_u_three)
      n_CSF3<-getFreq("CSF",n_CSF_u_three)
      n_D163<-getFreq("D16",n_D16_u_three)
      n_D53<-getFreq("D5", n_D5_u_three)
      n_vWA3<-getFreq("vWA",n_vWA_u_three)
      n_D33<-getFreq("D3",n_D3_u_three)
      n_D73<-getFreq("D7",n_D7_u_three)
      n_D83<-getFreq("D8",n_D8_u_three)
      n_D183<-getFreq("D18",n_D18_u_three)
      n_FGA3<-getFreq("FGA",n_FGA_u_three)
      n_TH013<-getFreq("TH01",n_TH01_u_three)
      n_D23<-getFreq("D2",n_D2_u_three)
      n_D193<-getFreq("D19",n_D19_u_three)
      n_D213<-getFreq("D21",n_D21_u_three)
      n_D133<-getFreq("D13",n_D13_u_three)
      
      
      rep_one_TPOX<-as.numeric(input$rep_one_TPOX)
      #print(rep_one_TPOX)
      rep_one_CSF<-as.numeric(input$rep_one_CSF)
      rep_one_D3<-as.numeric(input$rep_one_D3)
      rep_one_D16<-as.numeric(input$rep_one_D16)
      rep_one_D7<-as.numeric(input$rep_one_D7)
      rep_one_D8<-as.numeric(input$rep_one_D8)
      rep_one_D21<-as.numeric(input$rep_one_D21)
      rep_one_D18<-as.numeric(input$rep_one_D18)
      rep_one_FGA<-as.numeric(input$rep_one_FGA)
      rep_one_D5<-as.numeric(input$rep_one_D5)
      rep_one_D13<-as.numeric(input$rep_one_D13)
      rep_one_vWA<-as.numeric(input$rep_one_vWA)
      rep_one_TH01<-as.numeric(input$rep_one_TH01)
      rep_one_D2<-as.numeric(input$rep_one_D2)
      rep_one_D19<-as.numeric(input$rep_one_D19)
      
      
      rep_two_TPOX<-as.numeric(input$rep_two_TPOX)
      rep_two_CSF<-as.numeric(input$rep_two_CSF)
      rep_two_D3<-as.numeric(input$rep_two_D3)
      rep_two_D16<-as.numeric(input$rep_two_D16)
      rep_two_D7<-as.numeric(input$rep_two_D7)
      rep_two_D8<-as.numeric(input$rep_two_D8)
      rep_two_D21<-as.numeric(input$rep_two_D21)
      rep_two_D18<-as.numeric(input$rep_two_D18)
      rep_two_FGA<-as.numeric(input$rep_two_FGA)
      rep_two_D5<-as.numeric(input$rep_two_D5)
      rep_two_D13<-as.numeric(input$rep_two_D13)
      rep_two_vWA<-as.numeric(input$rep_two_vWA)
      rep_two_TH01<-as.numeric(input$rep_two_TH01)
      rep_two_D2<-as.numeric(input$rep_two_D2)
      rep_two_D19<-as.numeric(input$rep_two_D19)
      
      
      rep_three_TPOX<-as.numeric(input$rep_three_TPOX)
      rep_three_CSF<-as.numeric(input$rep_three_CSF)
      rep_three_D3<-as.numeric(input$rep_three_D3)
      rep_three_D16<-as.numeric(input$rep_three_D16)
      rep_three_D7<-as.numeric(input$rep_three_D7)
      rep_three_D8<-as.numeric(input$rep_three_D8)
      rep_three_D21<-as.numeric(input$rep_three_D21)
      rep_three_D18<-as.numeric(input$rep_three_D18)
      rep_three_FGA<-as.numeric(input$rep_three_FGA)
      rep_three_D5<-as.numeric(input$rep_three_D5)
      rep_three_D13<-as.numeric(input$rep_three_D13)
      rep_three_vWA<-as.numeric(input$rep_three_vWA)
      rep_three_TH01<-as.numeric(input$rep_three_TH01)
      rep_three_D2<-as.numeric(input$rep_three_D2)
      rep_three_D19<-as.numeric(input$rep_three_D19)
      
      D_U1<-rbind(getFreq("TPOX",d_TPOX_u_one), getFreq("CSF", d_CSF_u_one), getFreq("D3", d_D3_u_one), getFreq("D16", d_D16_u_one), getFreq("D7", d_D7_u_one)
                  , getFreq("D8", d_D8_u_one), getFreq("D21", d_D21_u_one), getFreq("D18", d_D18_u_one), getFreq("FGA", d_FGA_u_one), getFreq("D5", d_D5_u_one)
                  , getFreq("D13", d_D13_u_one), getFreq("vWA", d_vWA_u_one), getFreq("TH01", d_TH01_u_one), getFreq("D2", d_D2_u_one), getFreq("D19", d_D19_u_one))
      #print(D_U1)
      
      D_U2<-rbind(getFreq("TPOX",d_TPOX_u_two), getFreq("CSF", d_CSF_u_two), getFreq("D3", d_D3_u_two), getFreq("D16", d_D16_u_two), getFreq("D7", d_D7_u_two)
                  , getFreq("D8", d_D8_u_two), getFreq("D21", d_D21_u_two), getFreq("D18", d_D18_u_two), getFreq("FGA", d_FGA_u_two), getFreq("D5", d_D5_u_two)
                  , getFreq("D13", d_D13_u_two), getFreq("vWA", d_vWA_u_two), getFreq("TH01", d_TH01_u_two), getFreq("D2", d_D2_u_two), getFreq("D19", d_D19_u_two))
     #print(D_U2)
      
      D_U3<-rbind(getFreq("TPOX",d_TPOX_u_three), getFreq("CSF", d_CSF_u_three), getFreq("D3", d_D3_u_three), getFreq("D16", d_D16_u_three), getFreq("D7", d_D7_u_three)
                  , getFreq("D8", d_D8_u_three), getFreq("D21", d_D21_u_three), getFreq("D18", d_D18_u_three), getFreq("FGA", d_FGA_u_three), getFreq("D5", d_D5_u_three)
                  , getFreq("D13", d_D13_u_three), getFreq("vWA", d_vWA_u_three), getFreq("TH01", d_TH01_u_three), getFreq("D2", d_D2_u_three), getFreq("D19", d_D19_u_three))
      #print(D_U3)
      
      
      N_U1<-rbind(getFreq("TPOX",n_TPOX_u_one), getFreq("CSF", n_CSF_u_one), getFreq("D3", n_D3_u_one), getFreq("D16", n_D16_u_one), getFreq("D7", n_D7_u_one)
                  , getFreq("D8", n_D8_u_one), getFreq("D21", n_D21_u_one), getFreq("D18", n_D18_u_one), getFreq("FGA", n_FGA_u_one), getFreq("D5", n_D5_u_one)
                  , getFreq("D13", n_D13_u_one), getFreq("vWA", n_vWA_u_one), getFreq("TH01", n_TH01_u_one), getFreq("D2", n_D2_u_one), getFreq("D19", n_D19_u_one))
      #print(N_U1)
      
      N_U2<-rbind(getFreq("TPOX",n_TPOX_u_two), getFreq("CSF", n_CSF_u_two), getFreq("D3", n_D3_u_two), getFreq("D16", n_D16_u_two), getFreq("D7", n_D7_u_two)
                  , getFreq("D8", n_D8_u_two), getFreq("D21", n_D21_u_two), getFreq("D18", n_D18_u_two), getFreq("FGA", n_FGA_u_two), getFreq("D5", n_D5_u_two)
                  , getFreq("D13", n_D13_u_two), getFreq("vWA", n_vWA_u_two), getFreq("TH01", n_TH01_u_two), getFreq("D2", n_D2_u_two), getFreq("D19", n_D19_u_two))
      #print(N_U2)
      
      N_U3<-rbind(getFreq("TPOX",n_TPOX_u_three), getFreq("CSF", n_CSF_u_three), getFreq("D3", n_D3_u_three), getFreq("D16", n_D16_u_three), getFreq("D7", n_D7_u_three)
                  , getFreq("D8", n_D8_u_three), getFreq("D21", n_D21_u_three), getFreq("D18", n_D18_u_three), getFreq("FGA", n_FGA_u_three), getFreq("D5", n_D5_u_three)
                  , getFreq("D13", n_D13_u_three), getFreq("vWA", n_vWA_u_three), getFreq("TH01", n_TH01_u_three), getFreq("D2", n_D2_u_three), getFreq("D19", n_D19_u_three))
      #print(N_U3)
      
      TPOX_alleles<-unique(c(rep_one_TPOX, rep_two_TPOX, rep_three_TPOX))
      CSF_alleles<-unique(c(rep_one_CSF, rep_two_CSF, rep_three_CSF))
      D3_alleles<-unique(c(rep_one_D3, rep_two_D3, rep_three_D3))
      D16_alleles<-unique(c(rep_one_D16, rep_two_D16, rep_three_D16))
      D7_alleles<-unique(c(rep_one_D7, rep_two_D7, rep_three_D7))
      D8_alleles<-unique(c(rep_one_D8, rep_two_D8, rep_three_D8))
      D21_alleles<-unique(c(rep_one_D21, rep_two_D21, rep_three_D21))
      D18_alleles<-unique(c(rep_one_D18, rep_two_D18, rep_three_D18))
      FGA_alleles<-unique(c(rep_one_FGA, rep_two_FGA, rep_three_FGA))
      D5_alleles<-unique(c(rep_one_D5, rep_two_D5, rep_three_D5))
      D13_alleles<-unique(c(rep_one_D13, rep_two_D13, rep_three_D13))
      vWA_alleles<-unique(c(rep_one_vWA, rep_two_vWA, rep_three_vWA))
      TH01_alleles<-unique(c(rep_one_TH01, rep_two_TH01, rep_three_TH01))
      D2_alleles<-unique(c(rep_one_D2, rep_two_D2, rep_three_D2))
      D19_alleles<-unique(c(rep_one_D19, rep_two_D19, rep_three_D19))
      
      TPOX<-getFreq("TPOX", TPOX_alleles)
      
      #print("allele frequencis")
      #print(TPOX)
      #print(get_allele_freq(TPOX, "V3"))
      #TPOX_alleles_freq<-getFreq("TPOX", TPOX_alleles)
      #print(TPOX_alleles_freq)
      #temp<-get_allele_freq_2(TPOX_alleles_freq, "V3")
      #print("temp and other stuff")
      #print(temp)
      #allele_combos<-generate_allele_combos(temp)
      #print("allele combos")
      #print(allele_combos)
      #pg_table<-generate_freq_table(allele_combos)
      #all_the_combos<-generate_contributor_combos(allele_combos, number_contributors)
      #print("full freq table")
      #print(get_full_freq_table(pg_table, all_the_combos))
      CSF_alleles_freq<-getFreq("CSF", CSF_alleles)
      #print(TH01_alleles_freq)
      D3_alleles_freq<-getFreq("D3", D3_alleles)
      D16_alleles_freq<-getFreq("D16", D16_alleles)
      D7_alleles_freq<-getFreq("D7", D7_alleles)
      D8_alleles_freq<-getFreq("D8", D8_alleles)
      D21_alleles_freq<-getFreq("D21", D21_alleles)
      D18_alleles_freq<-getFreq("D18", D18_alleles)
      FGA_alleles_freq<-getFreq("FGA", FGA_alleles)
      D5_alleles_freq<-getFreq("D5", D5_alleles)
      D13_alleles_freq<-getFreq("D13", D13_alleles)
      vWA_alleles_freq<-getFreq("vWA", vWA_alleles)
      TH01_alleles_freq<-getFreq("TH01", TH01_alleles)
      D2_alleles_freq<-getFreq("D2", D2_alleles)
      D19_alleles_freq<-getFreq("D19", D19_alleles)
      
      Rep_1_TPOX<-getFreq("TPOX", rep_one_TPOX)
      #print("Rep 1 TPOX")
      #print(Rep_1_TPOX)
      #Rep_1_TPOX_dropout_table<-generate_rep_dropout_table(allele_combos, dropout_freq_table, Rep_1_TPOX)
      #print(get_full_rep_table(Rep_1_TPOX_dropout_table, all_the_combos))
      #print("full allele table")
      #print(get_full_allele_table(allele_combos, all_the_combos))
      Rep_1_CSF<-getFreq("CSF", rep_one_CSF)
      
      Rep_1_D3<-getFreq("D3", rep_one_D3)
      Rep_1_D16<-getFreq("D16", rep_one_D16)
      Rep_1_D7<-getFreq("D7", rep_one_D7)
      Rep_1_D8<-getFreq("D8", rep_one_D8)
      Rep_1_D21<-getFreq("D21", rep_one_D21)
      Rep_1_D18<-getFreq("D18", rep_one_D18)
      Rep_1_FGA<-getFreq("FGA", rep_one_FGA)
      Rep_1_D5<-getFreq("D5", rep_one_D5)
      Rep_1_D13<-getFreq("D13", rep_one_D13)
      Rep_1_vWA<-getFreq("vWA", rep_one_vWA)
      Rep_1_TH01<-getFreq("TH01", rep_one_TH01)
      Rep_1_D2<-getFreq("D2", rep_one_D2)
      Rep_1_D19<-getFreq("D19", rep_one_D19)
      
      
      Rep_2_TPOX<-getFreq("TPOX", rep_two_TPOX)
      Rep_2_CSF<-getFreq("CSF", rep_two_CSF)
      Rep_2_D3<-getFreq("D3", rep_two_D3)
      Rep_2_D16<-getFreq("D16", rep_two_D16)
      Rep_2_D7<-getFreq("D7", rep_two_D7)
      Rep_2_D8<-getFreq("D8", rep_two_D8)
      Rep_2_D21<-getFreq("D21", rep_two_D21)
      Rep_2_D18<-getFreq("D18", rep_two_D18)
      Rep_2_FGA<-getFreq("FGA", rep_two_FGA)
      Rep_2_D5<-getFreq("D5", rep_two_D5)
      Rep_2_D13<-getFreq("D13", rep_two_D13)
      Rep_2_vWA<-getFreq("vWA", rep_two_vWA)
      Rep_2_TH01<-getFreq("TH01", rep_two_TH01)
      Rep_2_D2<-getFreq("D2", rep_two_D2)
      Rep_2_D19<-getFreq("D19", rep_two_D19)
      
      
      Rep_3_TPOX<-getFreq("TPOX", rep_three_TPOX)
      Rep_3_CSF<-getFreq("CSF", rep_three_CSF)
      Rep_3_D3<-getFreq("D3", rep_three_D3)
      Rep_3_D16<-getFreq("D16", rep_three_D16)
      Rep_3_D7<-getFreq("D7", rep_three_D7)
      Rep_3_D8<-getFreq("D8", rep_three_D8)
      Rep_3_D21<-getFreq("D21", rep_three_D21)
      Rep_3_D18<-getFreq("D18", rep_three_D18)
      Rep_3_FGA<-getFreq("FGA", rep_three_FGA)
      Rep_3_D5<-getFreq("D5", rep_three_D5)
      Rep_3_D13<-getFreq("D13", rep_three_D13)
      Rep_3_vWA<-getFreq("vWA", rep_three_vWA)
      Rep_3_TH01<-getFreq("TH01", rep_three_TH01)
      Rep_3_D2<-getFreq("D2", rep_three_D2)
      Rep_3_D19<-getFreq("D19", rep_three_D19)
    
     #make_full_locus_table(D2_alleles_freq,Rep_1_D2, Rep_2_D2, Rep_3_D2, "D2", number_contributors, deducible_nondeducible, quantity, n_D21, n_D22,n_D23, d_D21, d_D22, d_D23, "V3")
     #make_full_locus_table(D13_alleles_freq,Rep_1_D13, Rep_2_D13, Rep_3_D13, "D13", number_contributors, deducible_nondeducible, quantity, n_D131, n_D132,n_D133, d_D131, d_D132, d_D133, "V3")
     #make_full_locus_table(CSF_alleles_freq,Rep_1_CSF, Rep_2_CSF, Rep_3_CSF, "CSF", number_contributors, deducible_nondeducible, quantity, n_CSF1, n_CSF2,n_CSF3, d_CSF1, d_CSF2, d_CSF3, "V3")
     make_full_locus_table(D7_alleles_freq,Rep_1_D7, Rep_2_D7, Rep_3_D7, "D7", number_contributors, deducible_nondeducible, quantity, n_D71, n_D72,n_D73, d_D71, d_D72, d_D73, "V3")
     #make_full_locus_table(D8_alleles_freq,Rep_1_D8, Rep_2_D8, Rep_3_D8, "D8", number_contributors, deducible_nondeducible, quantity, n_D81, n_D82,n_D83, d_D81, d_D82, d_D83, "V3")
     #make_full_locus_table(D21_alleles_freq,Rep_1_D21, Rep_2_D21, Rep_3_D21, "D21", number_contributors, deducible_nondeducible, quantity, n_D211, n_D212,n_D213, d_D211, d_D212, d_D213, "V3")
     #make_full_locus_table(D3_alleles_freq,Rep_1_D3, Rep_2_D3, Rep_3_D3, "D3", number_contributors, deducible_nondeducible, quantity, n_D31, n_D32,n_D33, d_D31, d_D32, d_D33, "V3")
     #make_full_locus_table(TH01_alleles_freq,Rep_1_TH01, Rep_2_TH01, Rep_3_TH01, "TH01", number_contributors, deducible_nondeducible, quantity, n_TH011, n_TH012,n_TH013, d_TH011, d_TH012, d_TH013, "V3")
     #make_full_locus_table(D16_alleles_freq,Rep_1_D16, Rep_2_D16, Rep_3_D16, "D16", number_contributors, deducible_nondeducible, quantity, n_D161, n_D162,n_D163, d_D161, d_D162, d_D163, "V3")
     #make_full_locus_table(D19_alleles_freq,Rep_1_D19, Rep_2_D19, Rep_3_D19, "D19", number_contributors, deducible_nondeducible, quantity, n_D191, n_D192,n_D193, d_D191, d_D192, d_D193, "V3")
     #make_full_locus_table(vWA_alleles_freq,Rep_1_vWA, Rep_2_vWA, Rep_3_vWA, "vWA", number_contributors, deducible_nondeducible, quantity, n_vWA1, n_vWA2,n_vWA3, d_vWA1, d_vWA2, d_vWA3, "V3")
     #make_full_locus_table(TPOX,Rep_1_TPOX, Rep_2_TPOX, Rep_3_TPOX, "TPOX", number_contributors, deducible_nondeducible, quantity, n_TPOX1, n_TPOX2,n_TPOX3, d_TPOX1, d_TPOX2, d_TPOX3, "V3")
     #make_full_locus_table(D18_alleles_freq,Rep_1_D18, Rep_2_D18, Rep_3_D18, "D18", number_contributors, deducible_nondeducible, quantity, n_D181, n_D182,n_D183, d_D181, d_D182, d_D183, "V3")
     #make_full_locus_table(D5_alleles_freq,Rep_1_D5, Rep_2_D5, Rep_3_D5, "D5", number_contributors, deducible_nondeducible, quantity, n_D51, n_D52,n_D53, d_D51, d_D52, d_D53, "V3")
     #make_full_locus_table(FGA_alleles_freq,Rep_1_FGA, Rep_2_FGA, Rep_3_FGA, "FGA", number_contributors, deducible_nondeducible, quantity, n_FGA1, n_FGA2,n_FGA3, d_FGA1, d_FGA2, d_FGA3, "V3")
     
     #print("Rep 1 D16")
     #print(Rep_1_vWA)
     #print("Rep 2 D16")
     #print(Rep_2_D16)
     #print("Rep 3 D16")
     #print(Rep_3_D16)
     #all_locuses<-rep(NA, 15)
     #all_locuses[1]<-make_full_locus_table(CSF_alleles_freq,Rep_1_CSF, Rep_2_CSF, Rep_3_CSF, "CSF", number_contributors, deducible_nondeducible, quantity, n_CSF1, n_CSF2,n_CSF3, d_CSF1, d_CSF2, d_CSF3, "V3")
     #all_locuses[2]<-make_full_locus_table(D3_alleles_freq,Rep_1_D3, Rep_2_D3, Rep_3_D3, "D3", number_contributors, deducible_nondeducible, quantity, n_D31, n_D32,n_D33, d_D31, d_D32, d_D33, "V3")
     #all_locuses[3]<-make_full_locus_table(D16_alleles_freq,Rep_1_D16, Rep_2_D16, Rep_3_D16, "D16", number_contributors, deducible_nondeducible, quantity, n_D161, n_D162,n_D163, d_D161, d_D162, d_D163, "V3")
     #all_locuses[4]<-make_full_locus_table(D7_alleles_freq,Rep_1_D7, Rep_2_D7, Rep_3_D7, "D7", number_contributors, deducible_nondeducible, quantity, n_D71, n_D72,n_D73, d_D71, d_D72, d_D73, "V3")
     #all_locuses[5]<-make_full_locus_table(D8_alleles_freq,Rep_1_D8, Rep_2_D8, Rep_3_D8, "D8", number_contributors, deducible_nondeducible, quantity, n_D81, n_D82,n_D83, d_D81, d_D82, d_D83, "V3")
     #all_locuses[6]<-make_full_locus_table(D21_alleles_freq,Rep_1_D21, Rep_2_D21, Rep_3_D21, "D21", number_contributors, deducible_nondeducible, quantity, n_D211, n_D212,n_D213, d_D211, d_D212, d_D213, "V3")
     #all_locuses[7]<-make_full_locus_table(D18_alleles_freq,Rep_1_D18, Rep_2_D18, Rep_3_D18, "D18", number_contributors, deducible_nondeducible, quantity, n_D181, n_D182,n_D183, d_D181, d_D182, d_D183, "V3")
     #all_locuses[8]<-make_full_locus_table(FGA_alleles_freq,Rep_1_FGA, Rep_2_FGA, Rep_3_FGA, "FGA", number_contributors, deducible_nondeducible, quantity, n_FGA1, n_FGA2,n_FGA3, d_FGA1, d_FGA2, d_FGA3, "V3")
     #all_locuses[9]<-make_full_locus_table(D5_alleles_freq,Rep_1_D5, Rep_2_D5, Rep_3_D5, "D5", number_contributors, deducible_nondeducible, quantity, n_D51, n_D52,n_D53, d_D51, d_D52, d_D53, "V3")
     #all_locuses[10]<-make_full_locus_table(D13_alleles_freq,Rep_1_D13, Rep_2_D13, Rep_3_D13, "D13", number_contributors, deducible_nondeducible, quantity, n_D131, n_D132,n_D133, d_D131, d_D132, d_D133, "V3")
     #all_locuses[11]<-make_full_locus_table(vWA_alleles_freq,Rep_1_vWA, Rep_2_vWA, Rep_3_vWA, "vWA", number_contributors, deducible_nondeducible, quantity, n_vWA1, n_vWA2,n_vWA3, d_vWA1, d_vWA2, d_vWA3, "V3")
     #all_locuses[12]<-make_full_locus_table(TH01_alleles_freq,Rep_1_TH01, Rep_2_TH01, Rep_3_TH01, "TH01", number_contributors, deducible_nondeducible, quantity, n_TH011, n_TH012,n_TH013, d_TH011, d_TH012, d_TH013, "V3")
     #all_locuses[13]<-make_full_locus_table(D2_alleles_freq,Rep_1_D2, Rep_2_D2, Rep_3_D2, "D2", number_contributors, deducible_nondeducible, quantity, n_D21, n_D22,n_D23, d_D21, d_D22, d_D23, "V3")
     #all_locuses[14]<-make_full_locus_table(D19_alleles_freq,Rep_1_D19, Rep_2_D19, Rep_3_D19, "D19", number_contributors, deducible_nondeducible, quantity, n_D191, n_D192,n_D193, d_D191, d_D192, d_D193, "V3")
     #all_locuses[15]<-make_full_locus_table(TPOX,Rep_1_TPOX, Rep_2_TPOX, Rep_3_TPOX, "TPOX", number_contributors, deducible_nondeducible, quantity, n_TPOX1, n_TPOX2,n_TPOX3, d_TPOX1, d_TPOX2, d_TPOX3, "V3")
     #print(all_locuses)
     #print(prod(as.numeric(all_locuses)))
      #TPOX_comb<-
     
      #z<- myFunction("D19", rep_three_D19)
      #print(z)
      
      })
   
  }
  
  shinyApp(ui = ui, server = server)
  
  
  #checkboxGroupInput(inputId = "numUsers", label = "Choose the number of users", choices = c("one" = "1", "two" = "2", "three" = "3", "four" = "4"),
                     #selected = "1", inline = FALSE)