

for i in doc-origin/*.html
do
   defname=`echo $i | sed "s/doc-origin\///"`
sed -e '1 s;^;<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"\
        "http://www.w3.org/TR/html4/loose.dtd">\
;' \
    -e '/Content-Type/ s;">;\; charset=UTF-8">;'\
    -e '/Content-Style-Type/ s;^;   <link rel="stylesheet" href="css/screen.css" type="text/css" media="screen, projection">\
   <link rel="stylesheet" href="css/print.css" type="text/css" media="print"> \
   <!--[if lt IE 8]>\
   <link rel="stylesheet" href="css/ie.css" type="text/css" media="screen, projection">\
   <![endif]-->\
   <link rel="stylesheet" href="css/texinfo.css" type="text/css" media="screen, projection"> \
   ;'\
    -e '/<body>/ s;$;\
    <div class="container">\
    <div class="span-24">;'\
    -e 's;<samp>;<code>;'\
    -e 's;</samp>;</code>;'\
    -e '/<\/body>/ s;$;\
    </div>\
    </div>;'\
    -e '/^<hr>$/ d' $i > $defname
done
