library(httr)
library(jsonlite)
library(dplyr)

# To use this script, you must be logged into COTSpotter in your browser.
# 1. Open your browser and log into https://cots.spot-lab.org
# 2. Open Developer Tools (F12) -> Application/Storage -> Cookies
# 3. Find the cookie named 'spotlab_token' and copy its value
# 4. Paste the value below:
# spotlab_token <- "YOUR_TOKEN_HERE"
spotlab_token <- "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJ1c2VybmFtZSI6InNhbSIsInByb2plY3RfaWQiOiJjb3RzIn0.sPawpJhWZb2Z4vWSTdH_Vwro64SWiEGPsB_riMgP3CM"
base_url <- "https://cots.spot-lab.org/api/external/v1/tracks"

# Initialize variables for pagination
all_tracks <- list()
current_offset <- 0
limit <- 100
keep_going <- TRUE

cat("Fetching COTS observations...\n")

while (keep_going) {
  # Construct the request
  response <- GET(
    url = base_url,
    query = list(
      aggregate = "all",
      limit = limit,
      offset = current_offset
    ),
    set_cookies(spotlab_token = spotlab_token)
  )

  # Check for errors
  if (status_code(response) == 401) {
    stop("Unauthorized. Please check if your spotlab_token is valid and not expired.")
  } else if (status_code(response) != 200) {
    stop(paste("API request failed with status", status_code(response)))
  }

  # Parse content
  content_text <- content(response, as = "text", encoding = "UTF-8")
  parsed_data <- fromJSON(content_text, flatten = TRUE)

  # Append data
  if (is.data.frame(parsed_data$data) && nrow(parsed_data$data) > 0) {
    all_tracks[[length(all_tracks) + 1]] <- parsed_data$data
  }

  # Check pagination
  next_offset <- parsed_data$pagination$next_offset

  cat(sprintf(
    "Fetched up to offset %d. Total available: %d\n",
    current_offset + nrow(parsed_data$data),
    parsed_data$pagination$total
  ))

  if (is.null(next_offset)) {
    keep_going <- FALSE
  } else {
    current_offset <- next_offset
    # Small delay to be polite to the server
    Sys.sleep(0.5)
  }
}

# Combine all pages into a single dataframe
cots_data <- bind_rows(all_tracks)

cat(sprintf("Successfully downloaded %d COTS observations.\n", nrow(cots_data)))

# Optionally save to a CSV
write.csv(cots_data, "reefscan_all_cots_observations.csv", row.names = FALSE)
