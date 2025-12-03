.onAttach <- function(libname, pkgname) {
  ver <- utils::packageVersion(pkgname)
  packageStartupMessage(paste0(
    "                          sss                     \n",
    "                 sss ssssssssssssss               \n",
    "            ssssssssssssssssssssssssssss          \n",
    "                 ssssssssssssssssssssssss         \n",
    "       sssssss ss   sssssssssssss s               \n",
    "         ssssssssssss   ss                        \n",
    "                 sss     ss sssssss               \n",
    "                      sssssssssssssssssssssssss   \n",
    "   ssssss s        sssssssssssssssssssssss  s     \n",
    "     ssssssssssss     ssssssssssssssssssss        \n",
    "                       ss    sssss                \n",
    "                     ssssssss                     \n",
    "                      ssssss                      \n",
    "                      ss   ss                     \n",
    "                     sss    ss                    \n",
    "                    sssssssssss                   \n",
    "                    ss       ss                   \n",
    "                    sssssssssss                   \n",
    "                    sss     ss                    \n",
    "                      ss  sss                     \n",
    "                         sss                      \n",
    "                       sss                        \n",
    "                     sss   sss                    \n",
    "                    sssssssssss                    \n",
    "                    ss       ss                   \n",
    "                    sssssssssss                   \n",
    "                    ss      sss       SequoiaR            \n",
    "                     ss    sss        Version: ", ver
  ))
}
