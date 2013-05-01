

for i in doc-origin/*.html
do
   defname=`echo $i | sed "s/doc-origin\///"`
sed -e '/style type=/ s;^;   <link rel="stylesheet" href="css/screen.css" type="text/css" media="screen, projection">\
   <link rel="stylesheet" href="css/print.css" type="text/css" media="print"> \
   <!--[if lt IE 8]>\
   <link rel="stylesheet" href="css/ie.css" type="text/css" media="screen, projection">\
   <![endif]-->\
   <link rel="stylesheet" href="css/texinfo.css" type="text/css" media="screen, projection"> \
   ;'\
    -e '/<body/ s;$;\
    <div class="container">\
    <div class="span-24">;'\
    -e 's;<samp>;<code>;'\
    -e 's;</samp>;</code>;'\
    -e '/<\/body>/ s;$;\
    </div>\
    </div>;'\
    -e '/^<hr>$/ d' $i > $defname
done

mv imogene.html index.html
