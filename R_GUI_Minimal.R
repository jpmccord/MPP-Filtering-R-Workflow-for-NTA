
filename <- tclvalue(tkgetOpenFile()) # Very simple, isn't it?
if (!nchar(filename)) {
  tkmessageBox(message = "No file was selected!")
} else {
  tkmessageBox(message = paste("The file selected was", filename))
}
  

tt <- tktoplevel()
cb <- tkcheckbutton(tt)
cbValue <- tclVar("0")
tkconfigure(cb, variable=cbValue)
tkgrid(tklabel(tt,text="I like R TclTk "),cb)
OnOK <- function()
{
  cbVal <- as.character(tclvalue(cbValue))
  tkdestroy(tt)
  if (cbVal=="1")
    tkmessageBox(message="So do I!")
  if (cbVal=="0")
    tkmessageBox(message="You forgot to check the box to say that you like R TclTk!",icon="warning")
}
OK.but <- tkbutton(tt,text="OK",command=OnOK)
tkgrid(OK.but)
tkfocus(tt)