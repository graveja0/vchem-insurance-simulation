if (!exists("manifest"))
    source(here::here("R/manifest.r"))

#https://asdfree.com/american-community-survey-acs.html

if (!file.exists("~/Desktop/acs-2012.rds")) {
    ####################
    # Household File
    ###################
    this_url_household <-
        "https://www2.census.gov/programs-surveys/acs/data/pums/2012/1-Year/csv_hus.zip"
    unzipped_files_household <- unzip("~/Downloads/csv_hus.zip" , exdir = tempdir())
    
    acs_csv_household <-
        grep('\\.csv$' , unzipped_files_household , value = TRUE)
    acs_df_household <- 
        acs_csv_household %>% map(~({
            data.table::fread(.x)
        })) %>%
        bind_rows()
        
    names(acs_df_household) <- tolower(names(acs_df_household))
    
    #################
    # Person File
    #################
    this_url_person <-
        "https://www2.census.gov/programs-surveys/acs/data/pums/2012/1-Year/csv_pus.zip"
    unzipped_files_person <- unzip("~/Downloads/csv_pus.zip" , exdir = tempdir())
    
    acs_csv_person <-
        grep('\\.csv$' , unzipped_files_person , value = TRUE)
    acs_df_person <- 
        acs_csv_person %>% map(~({
            data.table::fread(.x)
        })) %>% 
        bind_rows()
    names(acs_df_person) <- tolower(names(acs_df_person))
    
    acs_df_household[, 'rt'] <- NULL
    acs_df_person[, 'rt'] <- NULL
    acs_df <- merge(acs_df_household , acs_df_person)
    stopifnot(nrow(acs_df) == nrow(acs_df_person))
    acs_df[, 'one'] <- 1
    
    saveRDS(acs_df , file = "~/Desktop/acs-2012.rds" , compress = FALSE)
}

if (!file.exists("~/Desktop/acs-2022.rds")) {
    ####################
    # Household File
    ###################
    this_url_household <-
        "https://www2.census.gov/programs-surveys/acs/data/pums/2022/1-Year/csv_hus.zip"
    unzipped_files_household <- unzip("~/Downloads/csv_hus (1).zip" , exdir = tempdir())
    
    acs_csv_household <-
        grep('\\.csv$' , unzipped_files_household , value = TRUE)
    acs_df_household <- 
        acs_csv_household %>% map(~({
            data.table::fread(.x)
        })) %>%
        bind_rows()
    
    names(acs_df_household) <- tolower(names(acs_df_household))
    
    #################
    # Person File
    #################
    this_url_person <-
        "https://www2.census.gov/programs-surveys/acs/data/pums/2012/1-Year/csv_pus.zip"
    unzipped_files_person <- unzip("~/Downloads/csv_pus (1).zip" , exdir = tempdir())
    
    acs_csv_person <-
        grep('\\.csv$' , unzipped_files_person , value = TRUE)
    acs_df_person <- 
        acs_csv_person %>% map(~({
            data.table::fread(.x)
        })) %>% 
        bind_rows()
    names(acs_df_person) <- tolower(names(acs_df_person))
    
    acs_df_household[, 'rt'] <- NULL
    acs_df_person[, 'rt'] <- NULL
    acs_df <- merge(acs_df_household , acs_df_person)
    stopifnot(nrow(acs_df) == nrow(acs_df_person))
    acs_df[, 'one'] <- 1
    
    saveRDS(acs_df , file = "~/Desktop/acs-2022.rds" , compress = FALSE)
}

which_year = "2012"


acs_df <- readRDS(glue::glue("~/Desktop/acs-{which_year}.rds")) %>% 
    filter(agep<65) %>% 
    mutate(
        ins = if_else(hins1==1 | hins2==1 | hins3==1 | hins4==1 | hins5==1 | hins6==1, 1, 0),
        mcr = if_else(hins3==1, 1, 0),  # Medicare
        pub = if_else(hins4==1, 1, 0),  # Medicaid/means-tested (Other public)
        emp = if_else(hins1==1, 1, 0),  # Employer (maps to peg/poe)
        tri = if_else(hins5==1, 1, 0),  # TRICARE
        va = if_else(hins6==1, 1, 0),   # VA Health Care
        # Direct-purchase as non-group (maps to pne/png/pog/pri/prx)
        non_group = if_else(hins2==1, 1, 0),
        unin = 1 - ins
    ) %>%
    mutate(
        type = case_when(
            # 1) Insured & Medicare
            ins == 1 & mcr == 1 ~ 5,
            # 2) Insured & Other Public (not Medicare)
            ins == 1 & pub == 1 ~ 3,
            # 3) Insured & Employer/TRICARE/VA
            ins == 1 & (emp == 1 | tri == 1 | va == 1) ~ 1,
            # 4) Insured & Non-Group
            ins == 1 & non_group == 1 ~ 2,
            # 5) Uninsured
            unin == 1 ~ 4,
            TRUE ~ NA_real_
        )
    ) %>% 
    # Exclusion: No nonderly medicare!
    filter(type!=5)

acs_design <-
    svrepdesign(
        weight = ~pwgtp ,
        repweights = 'pwgtp[0-9]+' ,
        scale = 4 / 80 ,
        rscales = rep( 1 , 80 ) ,
        mse = TRUE ,
        type = 'JK1' ,
        data = acs_df
    )

acs_design <-
    update(
        
        acs_design ,
        
        state_name =
            factor(
                as.numeric( st ) ,
                levels = 
                    c(1L, 2L, 4L, 5L, 6L, 8L, 9L, 10L, 
                      11L, 12L, 13L, 15L, 16L, 17L, 18L, 
                      19L, 20L, 21L, 22L, 23L, 24L, 25L, 
                      26L, 27L, 28L, 29L, 30L, 31L, 32L, 
                      33L, 34L, 35L, 36L, 37L, 38L, 39L, 
                      40L, 41L, 42L, 44L, 45L, 46L, 47L, 
                      48L, 49L, 50L, 51L, 53L, 54L, 55L, 
                      56L, 72L) ,
                labels =
                    c("Alabama", "Alaska", "Arizona", "Arkansas", "California", 
                      "Colorado", "Connecticut", "Delaware", "District of Columbia", 
                      "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", 
                      "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", 
                      "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", 
                      "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", 
                      "New Mexico", "New York", "North Carolina", "North Dakota", "Ohio", 
                      "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina", 
                      "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", 
                      "Washington", "West Virginia", "Wisconsin", "Wyoming", "Puerto Rico")
            ) ,
        
        cit =
            factor( 
                cit , 
                levels = 1:5 , 
                labels = 
                    c( 
                        'born in the u.s.' ,
                        'born in the territories' ,
                        'born abroad to american parents' ,
                        'naturalized citizen' ,
                        'non-citizen'
                    )
            ) ,
        
        poverty_level = as.numeric( povpip ) ,
        
        married = as.numeric( mar %in% 1 ) ,
        
        sex = factor( sex , labels = c( 'male' , 'female' ) )
    )


acs_srvyr_design <- as_survey( acs_design )

# acs_srvyr_design %>% 
#     group_by(type) %>% summarise(mean = survey_mean())

res <- acs_srvyr_design %>%
    group_by(agep,type) %>% 
    summarise(pct = survey_mean())

res %>% 
    write_rds(here::here(glue::glue("_inputs/acs-calibration-targets/acs-calibration-targets_{which_year}.rds")))

res %>% 
    group_by(agep)  %>% 
    mutate(type = factor(type)) %>% 
    ggplot(aes(x = agep, y = pct,colour = type))  + geom_point() + geom_line() + geom_errorbar(aes(ymin = pct -1.96*pct_se,ymax=pct+1.96*pct_se))+
    ggthemes::theme_calc()


    read_rds(here::here("_inputs/acs-calibration-targets/acs-calibration-targets_2012.rds")) %>% 
    mutate(time = "pre") %>% 
    bind_rows(
        read_rds(here::here("_inputs/acs-calibration-targets/acs-calibration-targets_2022.rds")) %>% mutate(time = "post")
    ) %>% 
        mutate(type = factor(type)) %>% 
        ggplot(aes(x = agep, y = pct,colour = time))  + geom_point() + geom_line() + geom_errorbar(aes(ymin = pct -1.96*pct_se,ymax=pct+1.96*pct_se))+
        ggthemes::theme_calc() + facet_wrap(~type)

