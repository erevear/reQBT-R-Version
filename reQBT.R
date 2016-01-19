library(shiny)
frequencies<-read.csv("Allele_Frequencies.csv",header = FALSE, stringsAsFactors = FALSE)
drop_out_frequencies<-read.csv("Drop_Out_Rates.csv", header = FALSE, stringsAsFactors = FALSE)
colnames(drop_out_frequencies)<-c("Type", "Locus", 25, 50, 100, 150, 250, 500, 6.25, 12.5)
amounts<-c(NA, NA, 25,50,100,150,250,500, 6.25, 12.5)
#pn_known1<-c(12,14)
#pn_known2<-c()
#pn_known3<-c()
#dn_known1<-c()
#dn_known2<-c()
#dn_known3<-c()
#rep_one<-c(8,12,13,14,15)
#rep_two<-c(8,10,12,13,14,15)
#rep_three<-c(8,13,14)
#locus<-"D8"
#D_ND="D"
#quant<-60
#number_contributors<-3
#race<-"V3"
#rep_three<-c()
#o<-main_method(rep_one, rep_two, rep_three, quant, D_ND, number_contributor, locus, race, 
#                      pn_known1, pn_known2, pn_known3, dn_known1, dn_known2, dn_known3)
#o
get_present_alleles<-function(rep_one, rep_two, rep_three){
  return(unique(c(rep_one, rep_two, rep_three, 3.3)))
}



alleles

get_frequencies<-function(locus, alleles){
  y<-frequencies[frequencies$V1 == locus & frequencies$V2 %in% alleles, ]
 
  return(y)
}

connect_tables<-function(values, names, allele_vector){
  names(values)<-names
  vectorized <- values[as.character(allele_vector)]
  vectorized<-unname(vectorized)
  return(vectorized)
}

get_probabilities_table<-function(alleles, frequencies){
  alleles[,3]<-connect_tables(frequencies[,2], frequencies[,1], alleles[,1])
  alleles[,4]<-connect_tables(frequencies[,2], frequencies[,1], alleles[,2])
  alleles[,5] <- mapply(geno_freq, as.numeric(alleles[,1]), as.numeric(alleles[,2]), as.numeric(alleles[,3]), as.numeric(alleles[,4]))
  return(alleles)
}





geno_freq <- function(x,y,a,b){
  if(x==y){
    z<- a^2+0.03*a*(1-a)
  }
  else{
    z<- 2*a*b
  }
}

get_allele_freq_2<-function(alleles, race) {
  freq<-alleles[, race]
  temp<-as.numeric(alleles[,"V2"])
  the_names<-c(temp, 3.3)
  w<-1-sum(freq)
  freq[length(freq)+1]<-w
  names(freq)<-the_names
  return(freq)
}


#frequencies<-get_frequencies(locus, alleles)
#frequencies

get_allele_combos<-function(alleles){
  combos<-expand.grid(allele_1=alleles, allele_2=alleles)
  combos<-combos[combos[,2]>=combos[,1],] 
  
}

create_lookup_vector<-function(allele_combos){
  lookup_vector<-mapply(make_strings, allele_combos[,1], allele_combos[,2])
  return(lookup_vector)
}


make_strings<-function(allele_1, allele_2){
  string<-toString(c(allele_1, allele_2))
}

generate_all_combos<-function(allele_combos, num_contributors){
  all_user_combos<-expand.grid(replicate(num_contributors, allele_combos, simplify=FALSE))
  all_user_combos <- all_user_combos[,c(3, 2, 1)]
  return(all_user_combos)
}

get_dropout_range<-function(num_contributors, D_ND) {
  Het1_table<-get_dropout_types(num_contributors, D_ND, "HET1")
  Het2_table<-get_dropout_types(num_contributors, D_ND, "HET2")
  Hom1_table<-get_dropout_types(num_contributors, D_ND, "HOM1")
  full_table<-cbind(Hom1_table, Het1_table, Het2_table)
  return(full_table)

}

get_dropout_types<-function(num_contributors, D_ND, type){
  loc<-paste(type, num_contributors, D_ND, sep="-")
  table<-drop_out_frequencies[drop_out_frequencies$Type == loc,]
  return(table)
}

calculate_dropout<-function(dropout_rate_table, quant, locus) {
  temp_table<-dropout_rate_table[dropout_rate_table$Locus == locus,]
  possible_amounts<-c( 6.25, 12.5, 25,50,100,150,250,500)
  hom1<-temp_table[,1:10]
  het1<-temp_table[,11:20]
  het2<-temp_table[,21:30]

  if(quant%in%possible_amounts){
    Hom1_do_freq<-as.numeric(hom1[,toString(quant)])
    Het1_do_freq<-as.numeric(het1[,toString(quant)])
    Het2_do_freq<-as.numeric(het2[,toString(quant)])
  }else{
    
    possible_amounts[length(possible_amounts)+1]<-as.numeric(quant)
    possible_amounts<-sort(possible_amounts)
    
    Hom1_do_freq<-linear_interpolation(possible_amounts, quant, hom1)
    Het1_do_freq<-linear_interpolation(possible_amounts, quant, het1)
    Het2_do_freq<-linear_interpolation(possible_amounts, quant, het2)
  }
  Hom0_do_freq<- 1-Hom1_do_freq
  Het0_do_freq<- 1-(Het1_do_freq + Het2_do_freq)
  
  Hom_dropouts<-c(Hom1_do_freq,0, Hom0_do_freq)
  Het_dropouts<-c(Het2_do_freq, Het1_do_freq, Het0_do_freq)
  dropouts<-rbind(Hom_dropouts, Het_dropouts)

  return(dropouts) 
}

linear_interpolation<-function(amounts, quant, type_table){
  quant<-as.numeric(quant)
  quant_left<-amounts[which(amounts== as.numeric(quant)) - 1]
  quant_right<-amounts[which(amounts == as.numeric(quant)) + 1]
  type_quant_right<-as.numeric(type_table[,toString(quant_right)])
  type_quant_left<-as.numeric(type_table[,toString(quant_left)])

  type_slope<-(type_quant_right - type_quant_left)/(quant_right-quant_left)
  type_b_intercept<-type_quant_left - type_slope * quant_left
  type_do_freq<-type_slope * quant + type_b_intercept
  
  return(type_do_freq)
}

getFreq<-function(locus, alleles, race) {
  y<-frequencies[frequencies$V1 == locus & frequencies$V2 %in% alleles, ]
  y<-y[,2:3]
  w<-1-sum(y[,2])
  y[nrow(y)+1, 1]<-3.3
  y[nrow(y), 2]<-w
 
  return(y)
  
}

parse_table<-function(known1, known2, known3, full_table){
  allele_1<-c(known1[1], known2[1], known3[1])
  allele_2<-c(known1[2], known2[2], known3[2])
  print("here 4")
  pairs<-cbind(allele_1, allele_2)
  print(pairs)
  if(!is.null(pairs) & !all(is.na(pairs))){
    print("here 5")
    pairs_temp1<-pairs[,1]
    pairs_temp1<-pairs_temp1[!is.na(pairs_temp1)]
    pairs_temp2<-pairs[,2]
    pairs_temp2<-pairs_temp2[!is.na(pairs_temp2)]
    pairs<-cbind(pairs_temp1, pairs_temp2)
    #pairs[,1]<-pairs_temp1[!is.na(pairs_temp1)]
    #pairs[,2]<-pairs_temp2[!is.na(pairs_temp2)]
    print("pairs after")
    print(pairs)
    #pairs <- pairs[!is.na(pairs)]
    string_pairs<-mapply(make_strings, pairs[,1], pairs[,2])
  
    for(n in 1:length(string_pairs)){
      full_table<-full_table[full_table[,n]==string_pairs[n],]
    }
    full_table<-full_table[,(nrow(pairs) + 1):ncol(full_table)]
  }
  return(full_table)
}

check_known_rep<-function(known, allele_combos){
  if(length(known) >= 1){
    for(n in 1: length(known)){
      if(!(known[n]%in%allele_combos[,1]) & !(known[n]%in%allele_combos[,2])){
        known[n]<-3.3
        known<-known[order(known)]
      }
    }
  }
  return(known)
}
make_calculation<-function(full_table, pn_known1, pn_known2, pn_known3, dn_known1, dn_known2, dn_known3, number_contributors, allele_combos){
  print("here 3")
  pn_known1<-check_known_rep(pn_known1, allele_combos)
  pn_known2<-check_known_rep(pn_known2, allele_combos)
  pn_known3<-check_known_rep(pn_known3, allele_combos)
  dn_known1<-check_known_rep(dn_known1, allele_combos)
  dn_known2<-check_known_rep(dn_known2, allele_combos)
  dn_known3<-check_known_rep(dn_known3, allele_combos)
  print(pn_known3)
  print(dn_known2)
  numerator_table<-parse_table(pn_known1, pn_known2, pn_known3, full_table)
  numerator_table<-numerator_table[,(number_contributors + 1):ncol(numerator_table)]
  #numerator_table<-numerator_table[,(number_contributors + 1):ncol(numerator_table)]
  denominator_table<-parse_table(dn_known1, dn_known2, dn_known3, full_table)
  num<-apply(numerator_table, 1, prod)
  num<-sum(num)
  print(num)
  den<-apply(denominator_table[,(number_contributors + 1):ncol(denominator_table)], 1, prod)
  print(denominator_table)
  
  den<-sum(den)
  print(den)
  print(num/den)
  return(num/den)
}

modified_prod<-function(contributor_row, number_contributors){
  return(prod(contributor_row[,number_contributors:length(contributor_row)]))
}

apply_dropout<-function(alleles, dropout_rates, rep_info){
  alleles[,3]<-mapply(function(a, b){dropout<-0
                        
                               if(a == b){
                                 dropout<-dropout_rates[1, (1+(a%in%rep_info)+(b%in%rep_info))]
                               }else{
                                 dropout<-dropout_rates[2, (1+(a%in%rep_info)+(b%in%rep_info))]
                               }
                               return(dropout)}, alleles[,1], alleles[,2])
  return(alleles)
  
}

apply_dropin<-function(allele_combos, quant, rep_info){
  if(as.numeric(quant) <= 100){
    drop_in_frequencies<-c(PC0 = .96, PC1 = .035, PC2 = .005)
  }else{
    drop_in_frequencies<-c(PC0 = .975, PC1 = .02, PC2 = .005)
  }
  dropin<-apply(allele_combos, 1, function(allele_row){
      temp<-unique(as.numeric(unlist(strsplit(toString(as.matrix(allele_row)), split=","))))
      count<-length(which(temp%in%rep_info==TRUE))
      i<-1 + (length(rep_info) - count)
      ifelse(i <=3, drop_in<-drop_in_frequencies[i], drop_in<-drop_in_frequencies[3])
      return(drop_in)})
  return(dropin)
}

build_table<-function(rep_one, rep_two, rep_three, quant, D_ND, number_contributor, locus, race){
  alleles<-get_present_alleles(rep_one, rep_two, rep_three)
  print(alleles)
  allele_combos<-get_allele_combos(alleles)
  print(allele_combos)
  dropout_range<-get_dropout_range(number_contributors, D_ND)
  print(dropout_range)
  dropout_table<-calculate_dropout(dropout_range, quant, locus)
  print(dropout_table)
  allele_frequencies<-getFreq(locus, allele_combos[,1], race)
  print(allele_frequencies)
  probabilities<-get_probabilities_table(allele_combos, allele_frequencies)
  print(probabilities)
  rep_one_dropout<-apply_dropout(allele_combos, dropout_table, rep_one)
  rep_two_dropout<-apply_dropout(allele_combos, dropout_table, rep_two)
  rep_three_dropout<-apply_dropout(allele_combos, dropout_table, rep_three)
  print(rep_three_dropout)
  lookup_vector<-create_lookup_vector(allele_combos)
  all_combos<-generate_all_combos(lookup_vector, number_contributors)
  probabilities_cont1<-connect_tables(probabilities[,5], lookup_vector, all_combos[,1])
  print(probabilities_cont1)
  probabilities_cont2<-connect_tables(probabilities[,5], lookup_vector, all_combos[,2])
  probabilities_cont3<-connect_tables(probabilities[,5], lookup_vector, all_combos[,3])
  all_probabilities<-cbind(probabilities_cont1, probabilities_cont2, probabilities_cont3)
  u1_do<-get_all_dropout(all_combos[,1], rep_one_dropout, rep_two_dropout, rep_three_dropout, lookup_vector)
  u2_do<-get_all_dropout(all_combos[,2], rep_one_dropout, rep_two_dropout, rep_three_dropout, lookup_vector)
  u3_do<-get_all_dropout(all_combos[,3], rep_one_dropout, rep_two_dropout, rep_three_dropout, lookup_vector)
  
  rep_one_dropin<-apply_dropin(all_combos, quant, rep_one)
  rep_two_dropin<-apply_dropin(all_combos, quant, rep_two)
  rep_three_dropin<-apply_dropin(all_combos, quant, rep_three)
  return(cbind(all_combos, probabilities_cont1, probabilities_cont2, probabilities_cont3, u1_do, u2_do, u3_do, rep_one_dropin, rep_two_dropin, rep_three_dropin))
  
}

get_all_dropout<-function(alleles, do_rep1, do_rep2, do_rep3, lookup_vector){
  rep_one<-connect_tables(do_rep1[,3], lookup_vector, alleles)
  rep_two<-connect_tables(do_rep2[,3], lookup_vector, alleles)
  rep_three<-connect_tables(do_rep3[,3], lookup_vector, alleles)
  return(cbind(rep_one, rep_two, rep_three))
  
}

main_method<-function(rep_one, rep_two, rep_three, quant, D_ND, number_contributor, locus, race, 
                      pn_known1, pn_known2, pn_known3, dn_known1, dn_known2, dn_known3) {
  
  alleles<-get_present_alleles(rep_one, rep_two, rep_three)
  print(alleles)
  allele_combos<-get_allele_combos(alleles)
  r<-build_table(rep_one, rep_two, rep_three, quant, D_ND, number_contributor, locus, race)
  print("here")
  p<-make_calculation(r, pn_known1, pn_known2, pn_known3, dn_known1, dn_known2, dn_known3, number_contributors, allele_combos)
  return(p)

}

check_knowns<-function(knowns){
  if(length(knowns)==1){
    knowns[2]<-knowns[1]
  }
  return(knowns)
}

check_reps<-function(rep){
  if(length(rep) > 0) {
  
  
    if(rep[1]==1){
     rep<-c()
    }
    if(rep[1] == 2){
     
      rep[1]<-3.3
    }
  }
  return(rep)
}

#z<-get_allele_combos(alleles)
#z<-create_lookup_vector(z)
#z<-generate_all_combos(z, number_contributors)
#z
#x<-get_dropout_range(number_contributors, D_ND)
#x
#y<-calculate_dropout(x, quant, locus)
#y
#a<-getFreq(locus, z[,1], race)
#a
#b<-get_probabilities_table(z, a)
#b
#c<-apply_dropout(z, y, rep_one)
#c
#e<-create_lookup_vector(z)
#e
#d<-generate_all_combos(e, number_contributors)
#f<-connect_tables(b[,5], e, d[,1])
#g<-connect_tables(c[,3], e, d[, 1])
#h<-apply_dropin(d, quant, rep_one)
#h

ui<-shinyUI(fluidPage(
  
  mainPanel(
    tabsetPanel(
      tabPanel("Case Information", textInput(inputId = "Quantity", label = "Enter the quantity of DNA"), 
               numericInput(inputId = "numCont", label = "Select the number of contributors (1 - 4)", value = 1, min = 1, max = 4, step = NA),
               selectInput(inputId = "d_nd", label = "Is this profile deducible or non-deducible", choices = c("Deducible" = "D", "Non-Deducible" = "ND"), selected = "D", multiple = FALSE),
               actionButton(inputId = "submitButton", label = "Submit")),
      tabPanel("Denominator", tags$h3("Enter alleles of known contributors in the denominator"), tags$h4("K1"),
               checkboxGroupInput(inputId = "d_known_u_one_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D7", label = "D7",choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_vWA", label = "vWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_one_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
               
               
               
               
               
               
               tags$h4("K2"),
               
               checkboxGroupInput(inputId = "d_known_u_two_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_vWA", label = "vWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_two_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
               
               
               
               
               
               
               
               tags$h4("K3"),
               checkboxGroupInput(inputId = "d_known_u_three_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_vWA", label = "vWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "d_known_u_three_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE)
               
               
               
               
               
               
               
               #actionButton(inputId = "submitButton", label = "Submit")
               
      ),
      
      tabPanel("Numerator", tags$h3("Enter alleles of known contributors in the numerator"), tags$h4("K1"), 
               checkboxGroupInput(inputId = "n_known_u_one_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_one_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
               
               
               
               tags$h4("K2"),
               checkboxGroupInput(inputId = "n_known_u_two_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_two_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
               
               
               
               
               
               
               
               
               tags$h4("K3"),
               checkboxGroupInput(inputId = "n_known_three_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_three_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_three_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "n_known_u_three_FGA", label = "FGA",choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE)
               
               
               
               
               
               
               
               #actionButton(inputId = "submitButton", label = "Submit")
               
      ),
      
      tabPanel("Replicates", tags$h3("Enter replicate results"),tags$h4("Rep 1"), 
               checkboxGroupInput(inputId = "rep_one_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D21", label = "D21",choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_one_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
               
               
               
               
               tags$h4("Rep 2"),
               checkboxGroupInput(inputId = "rep_two_D8", label = "D8",choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_two_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE),
               
               
               
               tags$h4("Rep 3"),
               checkboxGroupInput(inputId = "rep_three_D8", label = "D8", choices = frequencies[frequencies$V1 == "D8", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D21", label = "D21", choices = frequencies[frequencies$V1 == "D21", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D7", label = "D7", choices = frequencies[frequencies$V1 == "D7", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_CSF", label = "CSF", choices = frequencies[frequencies$V1 == "CSF", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D3", label = "D3", choices = frequencies[frequencies$V1 == "D3", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_TH01", label = "TH01", choices = frequencies[frequencies$V1 == "TH01", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D13", label = "D13", choices = frequencies[frequencies$V1 == "D13", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D16", label = "D16", choices = frequencies[frequencies$V1 == "D16", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D2", label = "D2", choices = frequencies[frequencies$V1 == "D2", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D19", label = "D19", choices = frequencies[frequencies$V1 == "D19", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_vWA", label = "VWA", choices = frequencies[frequencies$V1 == "vWA", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_TPOX", label = "TPOX", choices = frequencies[frequencies$V1 == "TPOX", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D18", label = "D18", choices = frequencies[frequencies$V1 == "D18", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_D5", label = "D5", choices = frequencies[frequencies$V1 == "D5", "V2"], inline = TRUE),
               checkboxGroupInput(inputId = "rep_three_FGA", label = "FGA", choices = frequencies[frequencies$V1 == "FGA", "V2"], inline = TRUE)
               
               
               
               
               
               
               
               
               
      ),
      
      tabPanel("Results", textOutput("label_1"), tableOutput("suspect_summary"), tableOutput("suspect_loci"), tableOutput("false_positives"))
      
    ))
  

))


server<- function(input, output){
  observeEvent(input$submitButton, {
    quantity<-input$Quantity
    number_contributors<-input$numCont
    deducible_nondeducible<-input$d_nd
    
  
    d_TPOX_u_one<-check_knowns(as.numeric(input$d_known_u_one_TPOX))
    d_CSF_u_one<-check_knowns(as.numeric(input$d_known_u_one_CSF))
    d_D3_u_one<-check_knowns(as.numeric(input$d_known_u_one_D3))
    d_D16_u_one<-check_knowns(as.numeric(input$d_known_u_one_D16))
    d_D7_u_one<-check_knowns(as.numeric(input$d_known_u_one_D7))
    d_D8_u_one<-check_knowns(as.numeric(input$d_known_u_one_D8))
    d_D21_u_one<-check_knowns(as.numeric(input$d_known_u_one_D21))
    d_D18_u_one<-check_knowns(as.numeric(input$d_known_u_one_D18))
    d_FGA_u_one<-check_knowns(as.numeric(input$d_known_u_one_FGA))
    d_D5_u_one<-check_knowns(as.numeric(input$d_known_u_one_D5))
    d_D13_u_one<-check_knowns(as.numeric(input$d_known_u_one_D13))
    d_vWA_u_one<-check_knowns(as.numeric(input$d_known_u_one_vWA))
    d_TH01_u_one<-check_knowns(as.numeric(input$d_known_u_one_TH01))
    d_D2_u_one<-check_knowns(as.numeric(input$d_known_u_one_D2))
    d_D19_u_one<-check_knowns(as.numeric(input$d_known_u_one_D19))
    
    
    d_TPOX_u_two<-check_knowns(as.numeric(input$d_known_u_two_TPOX))
    d_CSF_u_two<-check_knowns(as.numeric(input$d_known_u_two_CSF))
    d_D3_u_two<-check_knowns(as.numeric(input$d_known_u_two_D3))
    d_D16_u_two<-check_knowns(as.numeric(input$d_known_u_two_D16))
    d_D7_u_two<-check_knowns(as.numeric(input$d_known_u_two_D7))
    d_D8_u_two<-check_knowns(as.numeric(input$d_known_u_two_D8))
    d_D21_u_two<-check_knowns(as.numeric(input$d_known_u_two_D21))
    d_D18_u_two<-check_knowns(as.numeric(input$d_known_u_two_D18))
    d_FGA_u_two<-check_knowns(as.numeric(input$d_known_u_two_FGA))
    d_D5_u_two<-check_knowns(as.numeric(input$d_known_u_two_D5))
    d_D13_u_two<-check_knowns(as.numeric(input$d_known_u_two_D13))
    d_vWA_u_two<-check_knowns(as.numeric(input$d_known_u_two_vWA))
    d_TH01_u_two<-check_knowns(as.numeric(input$d_known_u_two_TH01))
    d_D2_u_two<-check_knowns(as.numeric(input$d_known_u_two_D2))
    d_D19_u_two<-check_knowns(as.numeric(input$d_known_u_two_D19))
    
    
    d_TPOX_u_three<-check_knowns(as.numeric(input$d_known_u_three_TPOX))
    d_CSF_u_three<-check_knowns(as.numeric(input$d_known_u_three_CSF))
    d_D3_u_three<-check_knowns(as.numeric(input$d_known_u_three_D3))
    d_D16_u_three<-check_knowns(as.numeric(input$d_known_u_three_D16))
    d_D7_u_three<-check_knowns(as.numeric(input$d_known_u_three_D7))
    d_D8_u_three<-check_knowns(as.numeric(input$d_known_u_three_D8))
    d_D21_u_three<-check_knowns(as.numeric(input$d_known_u_three_D21))
    d_D18_u_three<-check_knowns(as.numeric(input$d_known_u_three_D18))
    d_FGA_u_three<-check_knowns(as.numeric(input$d_known_u_three_FGA))
    d_D5_u_three<-check_knowns(as.numeric(input$d_known_u_three_D5))
    d_D13_u_three<-check_knowns(as.numeric(input$d_known_u_three_D13))
    d_vWA_u_three<-check_knowns(as.numeric(input$d_known_u_three_vWA))
    d_TH01_u_three<-check_knowns(as.numeric(input$d_known_u_three_TH01))
    d_D2_u_three<-check_knowns(as.numeric(input$d_known_u_three_D2))
    d_D19_u_three<-check_knowns(as.numeric(input$d_known_u_three_D19))
    
  
    
    n_TPOX_u_one<-check_knowns(as.numeric(input$n_known_u_one_TPOX))
    n_CSF_u_one<-check_knowns(as.numeric(input$n_known_u_one_CSF))
    n_D3_u_one<-check_knowns(as.numeric(input$n_known_u_one_D3))
    n_D16_u_one<-check_knowns(as.numeric(input$n_known_u_one_D16))
    n_D7_u_one<-check_knowns(as.numeric(input$n_known_u_one_D7))
    n_D8_u_one<-check_knowns(as.numeric(input$n_known_u_one_D8))
    n_D21_u_one<-check_knowns(as.numeric(input$n_known_u_one_D21))
    n_D18_u_one<-check_knowns(as.numeric(input$n_known_u_one_D18))
    n_FGA_u_one<-check_knowns(as.numeric(input$n_known_u_one_FGA))
    n_D5_u_one<-check_knowns(as.numeric(input$n_known_u_one_D5))
    n_D13_u_one<-check_knowns(as.numeric(input$n_known_u_one_D13))
    n_vWA_u_one<-check_knowns(as.numeric(input$n_known_u_one_vWA))
    n_TH01_u_one<-check_knowns(as.numeric(input$n_known_u_one_TH01))
    n_D2_u_one<-check_knowns(as.numeric(input$n_known_u_one_D2))
    n_D19_u_one<-check_knowns(as.numeric(input$n_known_u_one_D19))
    
    
    n_TPOX_u_two<-check_knowns(as.numeric(input$n_known_u_two_TPOX))
    n_CSF_u_two<-check_knowns(as.numeric(input$n_known_u_two_CSF))
    n_D3_u_two<-check_knowns(as.numeric(input$n_known_u_two_D3))
    n_D16_u_two<-check_knowns(as.numeric(input$n_known_u_two_D16))
    n_D7_u_two<-check_knowns(as.numeric(input$n_known_u_two_D7))
    n_D8_u_two<-check_knowns(as.numeric(input$n_known_u_two_D8))
    n_D21_u_two<-check_knowns(as.numeric(input$n_known_u_two_D21))
    n_D18_u_two<-check_knowns(as.numeric(input$n_known_u_two_D18))
    n_FGA_u_two<-check_knowns(as.numeric(input$n_known_u_two_FGA))
    n_D5_u_two<-check_knowns(as.numeric(input$n_known_u_two_D5))
    n_D13_u_two<-check_knowns(as.numeric(input$n_known_u_two_D13))
    n_vWA_u_two<-check_knowns(as.numeric(input$n_known_u_two_vWA))
    n_TH01_u_two<-check_knowns(as.numeric(input$n_known_u_two_TH01))
    n_D2_u_two<-check_knowns(as.numeric(input$n_known_u_two_D2))
    n_D19_u_two<-check_knowns(as.numeric(input$n_known_u_two_D19))
    
    
    n_TPOX_u_three<-check_knowns(as.numeric(input$n_known_u_three_TPOX))
    n_CSF_u_three<-check_knowns(as.numeric(input$n_known_u_three_CSF))
    n_D3_u_three<-check_knowns(as.numeric(input$n_known_u_three_D3))
    n_D16_u_three<-check_knowns(as.numeric(input$n_known_u_three_D16))
    n_D7_u_three<-check_knowns(as.numeric(input$n_known_u_three_D7))
    n_D8_u_three<-check_knowns(as.numeric(input$n_known_u_three_D8))
    n_D21_u_three<-check_knowns(as.numeric(input$n_known_u_three_D21))
    n_D18_u_three<-check_knowns(as.numeric(input$n_known_u_three_D18))
    n_FGA_u_three<-check_knowns(as.numeric(input$n_known_u_three_FGA))
    n_D5_u_three<-check_knowns(as.numeric(input$n_known_u_three_D5))
    n_D13_u_three<-check_knowns(as.numeric(input$n_known_u_three_D13))
    n_vWA_u_three<-check_knowns(as.numeric(input$n_known_u_three_vWA))
    n_TH01_u_three<-check_knowns(as.numeric(input$n_known_u_three_TH01))
    n_D2_u_three<-check_knowns(as.numeric(input$n_known_u_three_D2))
    n_D19_u_three<-check_knowns(as.numeric(input$n_known_u_three_D19))
    
        
    rep_one_TPOX<-check_reps(as.numeric(input$rep_one_TPOX))
    rep_one_CSF<-check_reps(as.numeric(input$rep_one_CSF))
    rep_one_D3<-check_reps(as.numeric(input$rep_one_D3))
    rep_one_D16<-check_reps(as.numeric(input$rep_one_D16))
    rep_one_D7<-check_reps(as.numeric(input$rep_one_D7))
    rep_one_D8<-check_reps(as.numeric(input$rep_one_D8))
    rep_one_D21<-check_reps(as.numeric(input$rep_one_D21))
    rep_one_D18<-check_reps(as.numeric(input$rep_one_D18))
    rep_one_FGA<-check_reps(as.numeric(input$rep_one_FGA))
    rep_one_D5<-check_reps(as.numeric(input$rep_one_D5))
    rep_one_D13<-check_reps(as.numeric(input$rep_one_D13))
    rep_one_vWA<-check_reps(as.numeric(input$rep_one_vWA))
    rep_one_TH01<-check_reps(as.numeric(input$rep_one_TH01))
    rep_one_D2<-check_reps(as.numeric(input$rep_one_D2))
    rep_one_D19<-check_reps(as.numeric(input$rep_one_D19))
    
    
    rep_two_TPOX<-check_reps(as.numeric(input$rep_two_TPOX))
    rep_two_CSF<-check_reps(as.numeric(input$rep_two_CSF))
    rep_two_D3<-check_reps(as.numeric(input$rep_two_D3))
    rep_two_D16<-check_reps(as.numeric(input$rep_two_D16))
    rep_two_D7<-check_reps(as.numeric(input$rep_two_D7))
    rep_two_D8<-check_reps(as.numeric(input$rep_two_D8))
    rep_two_D21<-check_reps(as.numeric(input$rep_two_D21))
    rep_two_D18<-check_reps(as.numeric(input$rep_two_D18))
    rep_two_FGA<-check_reps(as.numeric(input$rep_two_FGA))
    rep_two_D5<-check_reps(as.numeric(input$rep_two_D5))
    rep_two_D13<-check_reps(as.numeric(input$rep_two_D13))
    rep_two_vWA<-check_reps(as.numeric(input$rep_two_vWA))
    rep_two_TH01<-check_reps(as.numeric(input$rep_two_TH01))
    rep_two_D2<-check_reps(as.numeric(input$rep_two_D2))
    rep_two_D19<-check_reps(as.numeric(input$rep_two_D19))
    
    
    rep_three_TPOX<-check_reps(as.numeric(input$rep_three_TPOX))
    rep_three_CSF<-check_reps(as.numeric(input$rep_three_CSF))
    rep_three_D3<-check_reps(as.numeric(input$rep_three_D3))
    rep_three_D16<-check_reps(as.numeric(input$rep_three_D16))
    rep_three_D7<-check_reps(as.numeric(input$rep_three_D7))
    rep_three_D8<-check_reps(as.numeric(input$rep_three_D8))
    rep_three_D21<-check_reps(as.numeric(input$rep_three_D21))
    rep_three_D18<-check_reps(as.numeric(input$rep_three_D18))
    rep_three_FGA<-check_reps(as.numeric(input$rep_three_FGA))
    rep_three_D5<-check_reps(as.numeric(input$rep_three_D5))
    rep_three_D13<-check_reps(as.numeric(input$rep_three_D13))
    rep_three_vWA<-check_reps(as.numeric(input$rep_three_vWA))
    rep_three_TH01<-check_reps(as.numeric(input$rep_three_TH01))
    rep_three_D2<-check_reps(as.numeric(input$rep_three_D2))
    rep_three_D19<-check_reps(as.numeric(input$rep_three_D19))
    
    
    #main_method<-function(rep_one, rep_two, rep_three, quant, D_ND, number_contributor, locus, race, 
     #                     pn_known1, pn_known2, pn_known3, dn_known1, dn_known2, dn_known3)
    
    results<-rep(NA, 4)
    races<-c("V3", "V4", "V5", "V6")
    #for(n in 1:4){
     race<-"V3"
      #race<-races[n]
      all_locuses<-rep(NA, 15)
      all_locuses[1]<-main_method(rep_one_D8, rep_two_D8, rep_three_D8, quantity, deducible_nondeducible, number_contributors,
                                  "D8", race, n_D8_u_one, n_D8_u_two, n_D8_u_three, d_D8_u_one, d_D8_u_two, d_D8_u_three)
      
      all_locuses[2]<-main_method(rep_one_D21, rep_two_D21, rep_three_D21, quantity, deducible_nondeducible, number_contributors,
                                  "D21", race, n_D21_u_one, n_D21_u_two, n_D21_u_three, d_D21_u_one, d_D21_u_two, d_D21_u_three  )
      
      all_locuses[3]<-main_method(rep_one_D7, rep_two_D7, rep_three_D7, quantity, deducible_nondeducible, number_contributors,
                                  "D7", race, n_D7_u_one, n_D7_u_two, n_D7_u_three, d_D7_u_one, d_D7_u_two, d_D7_u_three  )
      
      
      all_locuses[4]<-main_method(rep_one_CSF, rep_two_CSF, rep_three_CSF, quantity, deducible_nondeducible, number_contributors,
                                  "CSF", race, n_CSF_u_one, n_CSF_u_two, n_CSF_u_three, d_CSF_u_one, d_CSF_u_two, d_CSF_u_three  )
      
      all_locuses[5]<-main_method(rep_one_D3, rep_two_D3, rep_three_D3, quantity, deducible_nondeducible, number_contributors,
                                  "D3",race, n_D3_u_one, n_D3_u_two, n_D3_u_three, d_D3_u_one, d_D3_u_two, d_D3_u_three  )
      
      
      all_locuses[6]<-main_method(rep_one_TH01, rep_two_TH01, rep_three_TH01, quantity, deducible_nondeducible, number_contributors,
                                  "TH01", race, n_TH01_u_one, n_TH01_u_two, n_TH01_u_three, d_TH01_u_one, d_TH01_u_two, d_TH01_u_three  )
      
      
      all_locuses[7]<-main_method(rep_one_D13, rep_two_D13, rep_three_D13, quantity, deducible_nondeducible, number_contributors,
                                  "D13",race,  n_D13_u_one, n_D13_u_two, n_D13_u_three, d_D13_u_one, d_D13_u_two, d_D13_u_three  )
      
      all_locuses[8]<-main_method(rep_one_D16, rep_two_D16, rep_three_D16, quantity, deducible_nondeducible, number_contributors,
                                  "D16",race, n_D16_u_one, n_D16_u_two, n_D16_u_three, d_D16_u_one, d_D16_u_two, d_D16_u_three  )
      
      all_locuses[9]<-main_method(rep_one_D2, rep_two_D2, rep_three_D2, quantity, deducible_nondeducible, number_contributors,
                                  "D2", race, n_D2_u_one, n_D2_u_two, n_D2_u_three, d_D2_u_one, d_D2_u_two, d_D2_u_three  )
      
      all_locuses[10]<-main_method(rep_one_D19, rep_two_D19, rep_three_D19, quantity, deducible_nondeducible, number_contributors,
                                   "D19", race, n_D19_u_one, n_D19_u_two, n_D19_u_three, d_D19_u_one, d_D19_u_two, d_D19_u_three  )
      
      
      all_locuses[11]<-main_method(rep_one_vWA, rep_two_vWA, rep_three_vWA, quantity, deducible_nondeducible, number_contributors,
                                   "vWA", race, n_vWA_u_one, n_vWA_u_two, n_vWA_u_three, d_vWA_u_one, d_vWA_u_two, d_vWA_u_three  )
      
      
      all_locuses[12]<-main_method(rep_one_TPOX, rep_two_TPOX, rep_three_TPOX, quantity, deducible_nondeducible, number_contributors,
                                   "TPOX", race, n_TPOX_u_one, n_TPOX_u_two, n_TPOX_u_three, d_TPOX_u_one, d_TPOX_u_two, d_TPOX_u_three  )
      
      
      all_locuses[13]<-main_method(rep_one_D18, rep_two_D18, rep_three_D18, quantity, deducible_nondeducible, number_contributors,
                                   "D18", race, n_D18_u_one, n_D18_u_two, n_D18_u_three, d_D18_u_one, d_D18_u_two, d_D18_u_three  )
      
      
      all_locuses[14]<-main_method(rep_one_D5, rep_two_D5, rep_three_D5, quantity, deducible_nondeducible, number_contributors,
                                   "D5", race, n_D5_u_one, n_D5_u_two, n_D5_u_three, d_D5_u_one, d_D5_u_two, d_D5_u_three  )
      
      
      all_locuses[15]<-main_method(rep_one_FGA, rep_two_FGA, rep_three_FGA, quantity, deducible_nondeducible, number_contributors,
                                   "FGA", race, n_FGA_u_one, n_FGA_u_two, n_FGA_u_three, d_FGA_u_one, d_FGA_u_two, d_FGA_u_three  )
      
      
      print(all_locuses)
      print(prod(as.numeric(all_locuses)))
      #results[n]<-prod(as.numeric(all_locuses))
      
    #}
    #print("RESULTS")
    #print(results)
   
    
    
    
  
    
    
  
  })
  

    
    
}

shinyApp(ui = ui, server = server)


