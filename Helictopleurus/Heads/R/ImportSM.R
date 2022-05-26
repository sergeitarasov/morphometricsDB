# Function to Import stereoMorph ldks

importSM <- function(path, class = c('Opn', 'Out', 'Ldk'), panel = TRUE)
{
  require(magrittr)
  class <- match.arg(class)
  coo <- scale <- list()
  if (grepl('O', class)) {
    for (i in 1:length(list.files(path))) {
      parsed <- list.files(path, full.names = T)[i] %>% readLines
      coo.begin <- parsed %>% grep('<curves.pixel', .) + 1
      coo.end <- parsed %>% grep('</curves.pixel', .) - 1
      coo[[i]] <-
        parsed[coo.begin:coo.end] %>%
        gsub('^.*<.*>.*$', '', .) %>% #Removing separate curves
        sub('^\\t*', '', .) %>% #Removing tabs at the beginning and at the end of the strings
        strsplit('\t') %>%
        unlist %>%
        as.numeric %>%
        matrix(length(.)/2, 2, byrow = T)
      scale[i] <-
        parsed %>%
        extract(grep('<scaling[^.units]', .)) %>%
        sub('^.*>(.*)<.*$', '\\1', .) %>% #Leave only the number
        as.numeric
      names(coo)[i] <- list.files(path)[i]
    }
  } else {
    for (i in 1:length(list.files(path))) {
      parsed <- list.files(path, full.names = T)[i] %>% readLines
      coo.begin <- parsed %>% grep('<landmarks.pixel', .) + 1
      coo.end <- parsed %>% grep('</landmarks.pixel', .) - 1
      pre_ldk <-
        parsed[coo.begin:coo.end] %>%
        sub('^\\t*', '', .) %>% #Removing tabs at the beginning and at the end of the strings
        strsplit('\t')
      coo[[i]] <-
        pre_ldk %>%
        sapply(function(x) x[-1]) %>%
        apply(1, as.numeric) %>%
        set_rownames(pre_ldk %>% sapply(function(x) x[1])) %>%
        set_colnames(c('x', 'y')) %>%
        as.data.frame
      scale[i] <-
        parsed %>%
        extract(grep('<scaling[^.units]', .)) %>%
        sub('^.*>(.*)<.*$', '\\1', .) %>% #Leave only the number
        as.numeric
      names(coo)[i] <- list.files(path)[i]
    }
  }
  Coo <- match.arg(class) %>% do.call(list(coo, fac = data.frame(unlist(scale))))
  if (ncol(Coo$fac) != 0) names(Coo$fac) <- 'scale'
  if (panel) panel(Coo)
  return(Coo)
}
