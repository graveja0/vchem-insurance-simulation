if (!exists("manifest")) source(here::here("R/manifest.r"))

ins_types <- c("PRI","^INS.*X$","^MCR.*X$","^PUB.*X$","peg","poe","tri","vap","pne","png","pog","prx") %>% tolower()
# note "vap" converted to "va" in data processing steps!
final_types = c("pri","ins","mcr","pub","peg","poe","tri","va","pne","png","pog","prx")
final_types_ <- final_types[-which(final_types %in% c("poe","va","pne","prx"))] # used for pre data
type_lut <- c("type_1" = "emp_mil",
              "type_2" = "nongroup",
              "type_3" = "public",
              "type_4" = "uninsured",
              "type_5" = "medicare",
              "type_na" = "missing")

#-------------------------------------------------
# Download and save raw MEPS longitudinal files
#-------------------------------------------------

if (!file.exists(here("_data/meps/h245.rds"))) {
    # H245
    # Release date: November 2024
    # https://meps.ahrq.gov/mepsweb/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-245
    # Printing tipsThis file is a four-year longitudinal file derived from the respondents to the MEPS 
    # Panel 24 sample. The persons on this data set represent those who were in the MEPS population (U.S. 
    # civilian noninstitutionalized) for all or part of the 2019-2022 period. The file contains a 
    # longitudinal weight variable (LONGWT) and all variables from the 2019, 2020, 2021, and 2022 
    # full-year consolidated data files (HC-216, HC-224, HC-233, and HC-243, respectively). When 
    # applied to the persons who participated in 2019-2022, LONGWT enables users to make national 
    # estimates of person-level changes in selected variables (e.g., health insurance, health status, 
    # utilization, and expenditures). In addition, LONGWT can be used to develop cross-sectional type 
    # estimates for the four-year period and for each year individually based on only the Panel 24 sample. 
    meps_sas_import <- # Used for later years
        function( this_url ){
            
            this_tf <- tempfile()
            
            download.file( this_url , this_tf , mode = 'wb' )
            
            this_tbl <- read_sas( this_tf )
            
            this_df <- data.frame( this_tbl )
            
            names( this_df ) <- tolower( names( this_df ) )
            
            this_df
        }
    
    meps_post <-
        meps_sas_import( "https://meps.ahrq.gov/mepsweb/data_files/pufs/h245/h245v9.zip" )
    meps_post %>% write_rds(here("_data/meps/h245.rds"))
} else meps_post <- read_rds(here("_data/meps/h245.rds"))

if (!file.exists(here("_data/meps/h164.rds"))) {    
    # This file is a two-year longitudinal file derived from the respondents 
    # to the MEPS Panel 17 sample. The persons on this data set represent those 
    # who were in the MEPS population (U.S. civilian noninstitutionalized) for 
    # all or part of the 2012-2013 period. The file contains a longitudinal weight 
    # variable (LONGWT) and all variables from the 2012 and 2013 full-year consolidated 
    # data files (HC-155 and HC-163, respectively). The weight variable (LONGWT), when 
    # applied to the persons who participated in both 2012 and 2013, will enable the user to 
    # make national estimates of person-level changes in selected variables (e.g., health 
    # insurance, health status, utilization and expenditures). In addition, LONGWT can 
    # be used to develop cross-sectional type estimates for the two-year period and for
    # each year individually based on only the Panel 17 sample
    meps_sas_import_pre <- # Used for earlier years
        function( this_url ){
            
            this_tf <- tempfile()
            
            download.file( this_url , this_tf , mode = 'wb' )
            
            this_tbl <- read_xpt( this_tf )
            
            this_df <- data.frame( this_tbl )
            
            names( this_df ) <- tolower( names( this_df ) )
            
            this_df
        }
    meps_pre <- 
        meps_sas_import_pre("https://meps.ahrq.gov/mepsweb/data_files/pufs/h164ssp.zip")
    mepspre %>% write_rds(here("_data/meps/h164.rds"))
    
} else meps_pre <- read_rds(here("_data/meps/h164.rds"))

# Source: http://asdfree.com/medical-expenditure-panel-survey-meps.html
# ins_types <- c("GVA","GVB","GVC","HPD","HPE","HPN","HPO","HPR","HPX","IHS","^INS.*X$","^MCD.*X$","^MCR.*X$","PDK","PEG","PNE","PNG","POE","POG","PRI",
#                "PRX","^PUB.*X$","^TRI.*X$","^VAP.*(Y1|Y2|Y3|Y4)$") %>% tolower()
# final_types <- c("gva", "gvb", "gvc", "hpd", "hpe", "hpn", "hpo", "hpr", "hpx", "ihs", "ins", "mcd", "mcr", "pdk", "peg", "pne", "png", "poe", "pog", "pri", "prx", "pub", "tri", "vap")
# GVA{month}{year}: COV BY OTHER PUBLIC COVERAGE
# GVB{month}{year}: COV BY OTHER PUBLIC HMO
# GVC{month}{year}: COV BY OTHER PUBLIC PAYS PREM 
# HPD{month}{year}: PHOLDER OF PRIV INS (SRC UNKNWN) 
# HPE{month}{year}: PHOLDER OF EMPL UNION INS
# HPN{month}{year}: PHOLDER OF NONGROUP INS
# HPO{month}{year}: PHOLDER OF OTHER GROUP INS
# HPR{month}{year}: PHOLDER OF PRIVATE INSURANCE
# HPX{month}{year}: PHOLDER OF PRIV INS THRU EXCH
# IHS{month}{year}: COV BY INDIAN HEALTH SERVICE
# INS{month}{year}: COVR BY HOSP/MED INS 
# MCD{month}{year}{X}: COV BY MEDICAID OR SCHIP 
# MCR{month}{year}{X}: COV BY MEDICARE
# PDK{month}{year}: COVR BY PRIV INS (SOURCE UNKNWN) 
# PEG{month}{year}: COVERED BY EMPL UNION INS 
# PNE{month}{year}: COV BY NON-ESI,PHLDR OUTSIDE RU
# PNG{month}{year}: COVERED BY NONGROUP INS
# POE{month}{year}: COV BY ESI, PHOLDER OUTSIDE RU
# POG{month}{year}: COVERED BY OTHER GROUP INS
# PRI{month}{year}: COVERED BY PRIVATE INS
# PRX{month}{year}: COV BY PRIV INS THROUGH EXCHNG
# PUB{month}{year}{X}: COVR BY ANY PUBLIC INS 
# TRI{month}{year}{X}: COVERED BY TRICARE/CHAMPVA
# VAP{month}{year}{X}: COVERED BY VA

.x = ins_types[1]; .x

post_lut <- c("y1" = "2019",
              "y2" = "2020",
              "y3" = "2021",
              "y4" = "2022")
pre_lut <- c("y1" = "2012",
              "y2" = "2013",
              "y3" = "2014",
              "y4" = "2015") 
cal_lut <- c("ja" = "_01",
             "fe" = "_02",
             "ma" = "_03",
             "ap" = "_04",
             "my" = "_05",
             "ju" = "_06",
             "jl" = "_07",
             "au" = "_08",
             "se" = "_09",
             "oc" = "_10",
             "no" = "_11",
             "de" = "_12")
df_post <-
    ins_types %>% map( ~ ({
        cat(paste0(.x,"\n"))
        meps_post[, c("dupersid", grep(.x, colnames(meps_post), value = TRUE))] %>%
            as_tibble() %>%
            select(dupersid, matches("y[1-4]$"),matches("y[1-4]x$")) %>% {
                df <- .
                names(post_lut) %>% map( ~ (df <<- df %>% rename_all(function(x)
                    gsub(
                        .x, paste0(".", post_lut[.x]), x
                    ))))  
                
            } %>%
            
            pluck(length(post_lut)) %>%
            rename_all(function(x) gsub("x$","",x)) %>% 
            rename_all(function(x) gsub("^vap","vaa",x)) %>% 
            select(-contains("ev.")) %>% {
                df <- .
                names(cal_lut) %>% map( ~ (df <<- df %>% rename_all(function(x)
                    gsub(
                        .x, cal_lut[.x], x
                    ))))
            } %>%
            pluck((length(cal_lut))) %>%
            select(
                dupersid,
                contains("01."),
                contains("02."),
                contains("03."),
                contains("04."),
                contains("05."),
                contains("06."),
                contains("07."),
                contains("08."),
                contains("09."),
                contains("10."),
                contains("11."),
                contains("12.")
            ) %>%
            select(-contains(".p.")) %>%
            gather(tmp, value, -dupersid) %>%
            arrange(dupersid) %>%
            as_tibble() %>%
            mutate(
                month = as.numeric(substr(tmp, 5, 6)),
                year = as.numeric(substr(tmp, 8, 11)),
                type = substr(tmp, 1, 3)
            ) %>%
            mutate(type=ifelse(type=="vaa","va",type)) %>% 
            select(dupersid, type, month, year, value)  %>% 
            mutate(type = gsub("_","",type)) %>% 
            spread(type, value)
    }))

df_ins_post <- df_post[[1]]

tmp <- 2:length(df_post) %>% map(~({
    df_ins_post <<- df_ins_post %>% left_join(df_post[[.x]],c("dupersid","month","year"))
}))


df_pre <-
    ins_types[-which(ins_types %in% c("poe","vap","pne","prx"))] %>% map( ~ ({
        cat(paste0(.x,"\n"))
        meps_pre[, c("dupersid", grep(.x, colnames(meps_pre), value = TRUE))] %>%
            as_tibble() %>%
            select(dupersid, matches("y[1-4]$"),matches("y[1-4]x$")) %>% {
                df <- .
                names(pre_lut) %>% map( ~ (df <<- df %>% rename_all(function(x)
                    gsub(
                        .x, paste0(".", pre_lut[.x]), x
                    ))))  
                
            } %>%
            
            pluck(length(post_lut)) %>%
            rename_all(function(x) gsub("x$","",x)) %>% 
            rename_all(function(x) gsub("^vap","vaa",x)) %>% 
            select(-contains("ev.")) %>% {
                df <- .
                names(cal_lut) %>% map( ~ (df <<- df %>% rename_all(function(x)
                    gsub(
                        .x, cal_lut[.x], x
                    ))))
            } %>%
            pluck((length(cal_lut))) %>%
            select(
                dupersid,
                contains("01."),
                contains("02."),
                contains("03."),
                contains("04."),
                contains("05."),
                contains("06."),
                contains("07."),
                contains("08."),
                contains("09."),
                contains("10."),
                contains("11."),
                contains("12.")
            ) %>%
            select(-contains(".p.")) %>%
            gather(tmp, value, -dupersid) %>%
            arrange(dupersid) %>%
            as_tibble() %>%
            mutate(
                month = as.numeric(substr(tmp, 5, 6)),
                year = as.numeric(substr(tmp, 8, 11)),
                type = substr(tmp, 1, 3)
            ) %>%
            mutate(type=ifelse(type=="vaa","va",type)) %>% 
            select(dupersid, type, month, year, value)  %>% 
            mutate(type = gsub("_","",type)) %>% 
            spread(type, value)
    }))


df_ins_pre <- df_pre[[1]]

tmp <- 2:length(df_pre) %>% map(~({
    df_ins_pre <<- df_ins_pre %>% left_join(df_pre[[.x]],c("dupersid","month","year"))
}))

df_ins_post <- 
    df_ins_post %>% 
    mutate_at(vars(all_of(final_types)), function(x) ifelse(x==-1,NA,ifelse(x==2,0,ifelse(x==1,1,-99))))  %>% 
    mutate(unin = 1 - ins) %>% 
    mutate(type = case_when(
        ins == 1 & mcr == 1 ~ 5, # Medicare
        ins ==1 & pub == 1 ~ 3, # Public
        ins == 1 & (peg == 1 | poe == 1 | tri == 1 | va == 1) ~ 1, # Employer/TRICARE/VA
        ins == 1 & (pne == 1 | png == 1 | pog == 1 | pri == 1 | prx == 1) ~ 2, # Non-Group
        unin ==1 ~ 4, # Uninsured
        .default = NA
    )) 




df_baseline_post <- 
    meps_post %>% 
    filter(age1x!=-1) %>% 
    select(dupersid,dobmm,dobyy, starts_with("agey"),racev1x,racev2x,povlevy1,longwt)  %>% 
    as_tibble()

df_final_post <- 
    df_baseline_post %>% 
    inner_join(df_ins_post,"dupersid") %>% 
    mutate(age = year - dobyy + -1 * as.integer(month<dobmm)) %>% 
    select(-starts_with("agey"),-starts_with("dob")) %>% 
    arrange(dupersid, year, month) %>% 
    mutate(month = as.Date(str_c(year,str_pad(month,width=2,pad="0"),"01",sep="-")))

df_ins_pre <- 
    df_ins_pre %>% 
    mutate_at(vars(all_of(final_types_)), function(x) ifelse(x==-1,NA,ifelse(x==2,0,ifelse(x==1,1,-99))))  %>%
    mutate(unin = 1 - ins) %>% 
    mutate(type = case_when(
        ins == 1 & mcr == 1 ~ 5,
        ins ==1 & pub == 1 ~ 3,
        ins == 1 & (peg == 1  | tri == 1 ) ~ 1,
        ins == 1 & ( png == 1 | pog == 1 | pri == 1 ) ~ 2,
        
        unin == 1 ~ 4, 
        .default = NA
    )) 

df_baseline_pre <- 
    meps_pre %>% 
    filter(age1x!=-1) %>% 
    select(dupersid,dobmm,dobyy, starts_with("agey"),racev1x,racev2x, povlevy1, longwt) 

df_final_pre <-
    df_baseline_pre %>% 
    inner_join(df_ins_pre,"dupersid") %>% 
    mutate(age = year - dobyy + -1 * as.integer(month<dobmm)) %>% 
    select(-starts_with("agey"),-starts_with("dob")) %>% 
    arrange(dupersid, year, month) %>% 
    mutate(month = as.Date(str_c(year,str_pad(month,width=2,pad="0"),"01",sep="-")))


p_pre <- df_final_pre %>% 
    fastDummies::dummy_cols("type") %>% 
    group_by(month) %>% 
    summarise_at(vars(starts_with("type_")),~weighted.mean(.,w=longwt,na.rm=TRUE)) %>% 
    gather(type,value,-month)  %>% 
    {
        df <- .
        names(type_lut) %>% map(~(df <<- df %>% mutate(type = gsub(.x,type_lut[.x],type))))
    } %>% 
    pluck(length(type_lut)) %>% 
    ggplot(aes(x = month, y = value, colour = type)) + 
    geom_point() + 
    ggthemes::theme_clean() +
    scale_y_continuous(limits =c(0,1))

p_post <- df_final_post %>% 
    fastDummies::dummy_cols("type") %>% 
    group_by(month) %>% 
    summarise_at(vars(starts_with("type_")),~weighted.mean(.,w=longwt,na.rm=TRUE)) %>% 
    gather(type,value,-month)  %>% 
    {
        df <- .
        names(type_lut) %>% map(~(df <<- df %>% mutate(type = gsub(.x,type_lut[.x],type))))
    } %>% 
    pluck(length(type_lut)) %>% 
    ggplot(aes(x = month, y = value, colour = type)) + 
    geom_point() + 
    ggthemes::theme_clean() +
    scale_y_continuous(limits =c(0,1))

