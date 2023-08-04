#!/usr/bin/env Rscript

library(httr)

zscore <- function(data, url = "http://127.0.0.1:3030", id = 1) {
  start.time <- Sys.time()
  r <- POST(url, body = list(jsonrpc = "2.0", method = "zscore", params = data, id = id), encode = "json")
  end.time <- Sys.time()
  print(paste0("Time: ", end.time - start.time))

  if (status_code(r) == 200) {
    return(unlist(content(r, "parsed")$result))
  } else {
    return(c())
  }
}