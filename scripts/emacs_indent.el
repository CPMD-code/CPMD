; bash function used to cleaup: function cleanup { file=$1; grep 'DO *[0-9]\+' $file >&/dev/null; if [[ $? != 0 ]]; then emacs -batch $file -l /u/teo/Software/cpmd/CPMD/scripts/emacs_indent.el -f emacs-format-function ; fi ; /sp/fd/teo/Development/semd/trunk/tools/prettifier/prettify.py $file ; }
(defun emacs-format-function ()
   (f90-mode)
   (indent-region (point-min) (point-max) nil)
   (untabify (point-min) (point-max))
   (f90-upcase-keywords)
   (save-buffer)
)
