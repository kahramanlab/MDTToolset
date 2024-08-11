#' Find Isoform Switching Events
#'
#' @param cancer_ori MDT files of cancer cohort
#' @param gtex_ori MDT files of normal cohort
#' @param exp_values_pcawg abundance dataframe of cancer samples
#' @param exp_values_gtex abundance dataframe of gtex samples
#' @param ensp_sequences ensp sequences of transcripts
#' @param cutoff cutoff
#' @param cutoff_rate enrichment value (for example MDT1 should be expressed >=2 than the next one)
#' @param cutoff_other_MDTs cutoff to decide percentage of MDTs found in normal samples
#' @param ensp_only transcripts having ENSP ids
#'
#' @return datatable
#' @importFrom dplyr group_by mutate select count n arrange filter distinct desc summarise pull
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @import data.table
#' @export

switch_calculator <- function(cancer_ori, gtex_ori, exp_values_pcawg, exp_values_gtex, ensp_sequences, cutoff, cutoff_rate, cutoff_other_MDTs, ensp_only) {
  cancer <- cancer_ori %>% dplyr::filter(.data$rate >= cutoff_rate & .data$ENST1 %in% ensp_only$Transcript.stable.ID) # enst 2 is not important in this case - AK tool does not look into it
  gtex <- gtex_ori %>% dplyr::filter(.data$rate >= cutoff_rate & .data$ENST1 %in% ensp_only$Transcript.stable.ID)
  dMDT <- data.frame(SampleID = character(),
                     ENSG = character(),
                     dMDT = character(),
                     ENST2_cancer = character(),
                     TPM1_cancer = numeric(),
                     TPM2_cancer = numeric(),
                     enrichment = numeric(),
                     p_value = numeric(),
                     relative_cancer_exp = numeric(),
                     relative_gtex_exp = numeric(),
                     MDT_GTEx = character())
  dMDT3 <- data.frame(SampleID = character(),
                      ENSG = character(),
                      dMDT = character(),
                      ENST2_cancer = character(),
                      TPM1_cancer = numeric(),
                      TPM2_cancer = numeric(),
                      enrichment = numeric(),
                      p_value = numeric(),
                      relative_cancer_exp = numeric(),
                      relative_gtex_exp = numeric(),
                      MDT_GTEx = character())

  sample_number_gtex <- length(unique(gtex$SampleID)) # number of GTEx sample
  # Calculate the frequency of each ENST1 in gtex
  gtex_frequency <- gtex %>%
    dplyr::group_by(ENST1) %>%
    dplyr::summarise(Count = n())

  # Calculate the threshold number of samples
  threshold <- sample_number_gtex * (cutoff / 100)
  # Filter the frequencies based on the threshold
  filtered_gtex_ENST1 <- gtex_frequency %>%
    dplyr::filter(Count <= threshold) %>%
    dplyr::pull(ENST1)

 # Select ENST1s from cancer that are in the filtered list or not present in gtex at all
  potential_cMDT_lists <- cancer %>%
    dplyr::filter(ENST1 %in% filtered_gtex_ENST1 | !ENST1 %in% gtex$ENST1) %>%
    dplyr::pull(ENST1) %>%
    unique()

  for (i in potential_cMDT_lists) {
    statistical_test2 <- NULL
    statistical_test <- NULL
    gtex_exp_values_of_MDT <- NULL
    cancer_exp_values_of_MDT <- NULL
    MDT_highest_num <- NULL
    sequence_of_MDT <- NULL
    MDT_in_sample <- i
    data_of_MDT <- cancer %>% dplyr::filter(.data$ENST1 == MDT_in_sample)
    gene_id_sample <- data_of_MDT[1, 2]
    sample_ids <- data_of_MDT[, 1]
    l <- NULL

    MDT_list_in_gtex <- gtex %>% dplyr::filter(.data$ENST1 != MDT_in_sample & .data$ENSG == gene_id_sample)
    unique_gtex_ensts <- unique(MDT_list_in_gtex$ENST1)

    sequence_of_MDT <- ensp_sequences %>%
      dplyr::filter(.data$ENST == MDT_in_sample) %>%
      dplyr::select("ENST_Seq")

    for (each_enst in unique_gtex_ensts) {
      sequence_of_enst <- ensp_sequences %>%
        dplyr::filter(.data$ENST == each_enst) %>%
        dplyr::select("ENST_Seq")
      if (identical(sequence_of_enst$ENST_Seq, sequence_of_MDT$ENST_Seq)) {
        MDT_list_in_gtex <- MDT_list_in_gtex %>% dplyr::filter(.data$ENST1 != each_enst)
      }
    }


    MDT_highest <- MDT_list_in_gtex %>%
      dplyr::group_by(.data$ENST1) %>%
      dplyr::count() %>%
      dplyr::arrange(desc(n))
    MDT_highest_num <- as.numeric(MDT_highest[1, 2])

    unique_gtex_ensts_after_remove_redundant <- unique(MDT_list_in_gtex$ENST1)

    percent_of_MDT_in_GTEx_samples <- MDT_highest_num / sample_number_gtex * 100

    if (!is.na(percent_of_MDT_in_GTEx_samples)) {
      if (percent_of_MDT_in_GTEx_samples >= cutoff_other_MDTs) {
        number_of_MDT_in_gtex <- length(gtex[gtex$ENST1 == MDT_in_sample, 1])
        percent_of_MDT_in_gtex <- number_of_MDT_in_gtex / sample_number_gtex * 100

        cancer_relative_exp_of_MDT <- (exp_values_pcawg %>% dplyr::filter(.data$ENST == MDT_in_sample))[, sample_ids] / colSums(as.data.frame(exp_values_pcawg[exp_values_pcawg$ENSG == gene_id_sample, sample_ids]))

        gtex_relative_exp_of_MDT <- as.numeric((exp_values_gtex %>% dplyr::filter(.data$ENST == MDT_in_sample))[, 3:ncol(exp_values_gtex)] / colSums(exp_values_gtex[exp_values_gtex$ENSG == gene_id_sample, 3:ncol(exp_values_gtex)]))

        median_gtex_relative_exp_of_MDT <- median(as.numeric(as.vector(gtex_relative_exp_of_MDT)), na.rm = TRUE)

        if (is.data.frame(cancer_relative_exp_of_MDT)) {
          # If it has a single column, extract it as a vector
          if (ncol(cancer_relative_exp_of_MDT) == 1) {
            cancer_values <- cancer_relative_exp_of_MDT[[1]]
          } else if (nrow(cancer_relative_exp_of_MDT) == 1) {
            # If it has a single row, extract it as a vector
            cancer_values <- unlist(cancer_relative_exp_of_MDT[1, ])
          } else {
            stop("the data should be a single-column or single-row data frame.")
          }
        } else {
          # If it's not a data frame, assume it's a single value
          cancer_values <- as.vector(cancer_relative_exp_of_MDT)
        }

        statistical_test <- sapply(cancer_values, function(i) BSDA::SIGN.test(gtex_relative_exp_of_MDT, alternative = "less", md = i))["p.value", ]

        data_of_MDT_w_statistics <- cbind(data_of_MDT, unlist(statistical_test), t(cancer_relative_exp_of_MDT), replicate(nrow(data_of_MDT), median_gtex_relative_exp_of_MDT))

        colnames(data_of_MDT_w_statistics) <- c("SampleID", "ENSG", "ENST1_cancer", "ENST2_cancer", "TPM1_cancer", "TPM2_cancer", "enrichment", "p_value", "relative_cancer_exp", "relative_gtex_exp")

        data_of_MDT_w_statistics_filtered <- data_of_MDT_w_statistics %>%
            dplyr::filter(.data$p_value <= 0.05) %>%
            dplyr::filter(.data$relative_cancer_exp > median_gtex_relative_exp_of_MDT)

        if (isTRUE(!is.na(data_of_MDT_w_statistics_filtered$SampleID[1]))) {
            l <- cbind(data_of_MDT_w_statistics_filtered, replicate(nrow(data_of_MDT_w_statistics_filtered), paste(unique_gtex_ensts_after_remove_redundant, collapse = ",")))

          colnames(l) <- c("SampleID", "ENSG", "dMDT", "ENST2_cancer", "TPM1_cancer", "TPM2_cancer", "enrichment", "p_value", "relative_cancer_exp", "relative_gtex_exp", "MDT_GTEx")
          dMDT <- data.table::rbindlist(list(dMDT, l))
          }
      }
    }
  }

  if(nrow(dMDT) != 0){

    colnames(dMDT) <- c("SampleID", "ENSG", "dMDT", "ENST2_cancer", "TPM1_cancer", "TPM2_cancer", "enrichment", "p_value", "relative_cancer_exp", "relative_gtex_exp", "MDT_GTEx")
    dMDT <- as.data.frame(dMDT)
    dMDT2 <- dMDT[order(dMDT$p_value), ]
    dMDT3 <- cbind(dMDT2, stats::p.adjust(dMDT2$p_value, method = "BH"))
    colnames(dMDT3) <- c("SampleID", "ENSG", "dMDT", "ENST2_cancer", "TPM1_cancer", "TPM2_cancer", "enrichment", "p_value", "relative_cancer_exp", "relative_gtex_exp", "MDT_GTEx" , "adj_p")

  } else {
      dMDT3 <- print('There is no MDT switching events between datasets')
    }
  return(dMDT3)
}
