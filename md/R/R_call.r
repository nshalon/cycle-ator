#!/usr/local/bin Rscript

options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <r file> <function> [param1=v11 v12 ... param2=v21 v22 ...]\n", script.name))
  q(status=1)
}

source.fn = args[1]
func = args[2]
params = args[3:length(args)]
params = paste(params, collapse=" ")

s = sapply(strsplit(params, '\\s+', perl=T), function(x) x[x != ""])
s = strsplit(s, '=', perl=T)
keys = sapply(s,function(x) { x[1] })
values = sapply(s,function(x) { x[2] })

param.list = list()
for (i in seq_along(keys)) {
  key = keys[i]
  value = values[[i]]

  # first try to convert to numerical
  options(warn=-1) # to disable warning 'NAs introduced by coercion'
  if (all(!is.na(as.numeric(value))))
    value = as.numeric(value)
  options(warn=1)

  # check if we have booleans
  if (all(is.element(value, c("T", "F"))))
    value = (value == "T")

  param.list[[key]] = value
}

# load source file
cat(sprintf("loading %s\n", source.fn))
suppressPackageStartupMessages(source(source.fn))
options(error=NULL)

tostr = function(x)
{
  if (length(x) == 1 && x == "NULL")
    return (NULL)
  if (typeof(x) == "character")
    x = paste("\"", x, "\"", sep="")
  if (length(x) > 1)
    x = paste("c(", paste(x, collapse=",", sep=""), ")", sep="")
  x
}

# print the function call, helpful for debug
fmt.list = lapply(param.list, tostr)
str = paste(names(fmt.list), fmt.list, sep="=", collapse=", ")
cat(sprintf(">> Calling:\n%s(%s)\n", func, str))

str2 = paste(names(fmt.list), fmt.list, sep="=", collapse="; ")
cat(sprintf(">> Parameters:\n%s\n", str2))

# finally we omit NULL values from list
param.list[sapply(param.list, function(x) length(x) == 1 && x == "NULL")] = NULL

# call the function
do.call(func, param.list)
