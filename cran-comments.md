## Submission

* Updated version 0.5.11

## Comment

* Forgot a test when an argument was missing

## Test environments

* local WIN-64x install, R 4.3.0

* win-builder devel 

* r-hub.io


## R CMD check results

* There were no ERRORs and no WARNINGs.

* Win-dev returned zero note:

* R-hub does not work anymore. I get this error
  Error in curl::curl_fetch_memory(url, handle = handle) : 
  SSL peer certificate or SSH remote key was not OK: [builder.r-hub.io] schannel: SEC_E_UNTRUSTED_ROOT (0x80090325) 

## Downstream dependencies

* This package does not rely on any other packages

* It requires Rtools4.3

