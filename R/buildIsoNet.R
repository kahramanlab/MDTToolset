#' Title
#'
#' @param dMDT_network_info result of switch_calculator
#' @param ensembl_data ensembl data
#' @param iso_network isoform interaction network
#' @param kallisto_counts kallisto count of the network to check if the interaction partners are expressed
#'
#' @return dataframe
#' @importFrom dplyr group_by mutate select count n arrange filter distinct desc everything
#' @importFrom tidyr separate_rows
#' @export

buildIsoNet <- function (dMDT_network_info, ensembl_data, iso_network, kallisto_counts)
{
  dMDT_int_Loss_all <- data.frame()
  for (each_sample in unique(dMDT_network_info$SampleID)) {
    sample_kallisto <- kallisto_counts[, c("ENST", "ENSG",
                                           each_sample)]
    sample_dMDTs <- dMDT_network_info[dMDT_network_info$SampleID ==
                                        each_sample, ]
    for (each_dMDT in unique(sample_dMDTs$dMDT)) {
      ensp_expressed_missed <- c()
      ensp_expressed_kept <- c()
      each_dMDT_in_sample <- sample_dMDTs[sample_dMDTs$dMDT ==
                                            each_dMDT, ]
      miss_ints <- each_dMDT_in_sample[each_dMDT_in_sample$dMDT ==
                                         each_dMDT, "MissInts"]
      kept_ints <- each_dMDT_in_sample[each_dMDT_in_sample$dMDT ==
                                         each_dMDT, "ExistInts"]
      miss_protein_ints <- regmatches(miss_ints, gregexpr("ENSP\\d+",
                                                          miss_ints))[[1]]
      kept_protein_ints <- regmatches(kept_ints, gregexpr("ENSP\\d+",
                                                          kept_ints))[[1]]
      if (length(miss_protein_ints) > 0) {
        for (each_protein in miss_protein_ints) {
          ENST_id_missed <- unique(ensembl_data[ensembl_data$Protein.stable.ID ==
                                                  each_protein, "Transcript.stable.ID"])
          ENST_id_expression_missed <- sample_kallisto[sample_kallisto$ENST ==
                                                         ENST_id_missed, 3]
          if (length(ENST_id_expression_missed) > 0 &&
              as.numeric(ENST_id_expression_missed) >=
              2) {
            ensp_expressed_missed <- append(ensp_expressed_missed,
                                            each_protein)
          }
          else {
            ensp_expressed_missed <- append(ensp_expressed_missed,
                                            "Not_expressed")
          }
        }
      }
      else {
        ensp_expressed_missed <- append(ensp_expressed_missed,
                                        "No missed int")
      }
      if (length(kept_protein_ints) > 0) {
        for (each_protein_kept in kept_protein_ints) {
          ENST_id_kept <- unique(ensembl_data[ensembl_data$Protein.stable.ID ==
                                                each_protein_kept, "Transcript.stable.ID"])
          ENST_id_expression_kept <- sample_kallisto[sample_kallisto$ENST ==
                                                       ENST_id_kept, 3]
          if (length(ENST_id_expression_kept) > 0 &&
              as.numeric(ENST_id_expression_kept) >= 2) {
            ensp_expressed_kept <- append(ensp_expressed_kept,
                                          each_protein_kept)
          }
          else {
            ensp_expressed_kept <- append(ensp_expressed_kept,
                                          "Not_expressed")
          }
        }
      }
      else {
        ensp_expressed_kept <- append(ensp_expressed_kept,
                                      "No kept int")
      }
      MDT_in_samples <- each_dMDT_in_sample %>% tidyr::separate_rows(MDT_GTEx,
                                                                     sep = ",")
      # MDT_in_sample <- unique(each_dMDT_in_sample[, "MDT_GTEx"])
      MDT_in_sample <- unique(unlist(strsplit(each_dMDT_in_sample$MDT_GTEx, ",")))
      for (eachMDT in MDT_in_sample) {
        MDT_ENSP <- unique(ensembl_data[ensembl_data$Transcript.stable.ID ==
                                          eachMDT, "Protein.stable.ID"])
        Canonical_STRING_ENSP <- each_dMDT_in_sample[,
                                                     "STRINGensp"]
        if (length(MDT_ENSP) > 0) {
          if (MDT_ENSP != Canonical_STRING_ENSP) {
            MDT_int_losts <- iso_network[iso_network$ENST ==
                                           eachMDT, "MissInts"]
            MDT_int_losts_ensps <- regmatches(miss_ints,
                                              gregexpr("ENSP\\d+", MDT_int_losts))[[1]]
            int_losts_in_sample <- ensp_expressed_missed[!ensp_expressed_missed %in%
                                                           MDT_int_losts_ensps]
            int_losts_in_mdt <- ensp_expressed_missed[ensp_expressed_missed %in%
                                                        MDT_int_losts_ensps]
            each_dMDT_in_sample$ensp_expressed_missed <- unique(paste(ensp_expressed_missed,
                                                                      collapse = ","))
            each_dMDT_in_sample$ensp_expressed_kept <- unique(paste(ensp_expressed_kept,
                                                                    collapse = ","))
            each_dMDT_in_sample$int_losts_in_mdt <- unique(paste(int_losts_in_mdt,
                                                                 collapse = ","))
          }
          else if (MDT_ENSP == Canonical_STRING_ENSP) {
            each_dMDT_in_sample$ensp_expressed_missed <- unique(paste(ensp_expressed_missed,
                                                                      collapse = ","))
            each_dMDT_in_sample$ensp_expressed_kept <- unique(paste(ensp_expressed_kept,
                                                                    collapse = ","))
            each_dMDT_in_sample$int_losts_in_mdt <- "Canonical"
          }
        }
        else {
          each_dMDT_in_sample$ensp_expressed_missed <- unique(paste(ensp_expressed_missed,
                                                                    collapse = ","))
          each_dMDT_in_sample$ensp_expressed_kept <- unique(paste(ensp_expressed_kept,
                                                                  collapse = ","))
          each_dMDT_in_sample$int_losts_in_mdt <- "MDT_is_not_found"
        }
        dMDT_int_Loss_all <- rbind(dMDT_int_Loss_all,
                                   each_dMDT_in_sample)
      }
    }
  }
  return(dMDT_int_Loss_all)
}
