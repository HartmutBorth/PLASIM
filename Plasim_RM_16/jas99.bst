% Attempted hack to create BibTeX `jas' style for BibTeX version 0.99c
% Journal of Atmospheric Science, Journal of Applied Meteorology
% Monthly Weather Review
% by Mike Moran, Dept. of Atmospheric Science, Colorado State University
% based on `newapa' style file.  Changes are documented.  March 9, 1990
% Submitted 20-JUL-1993 by Piotr J. Flatau <pflatau@macao.ucsd.edu>
% e-mail for Mike Moran: <mmoran@cid.aes.doe.ca>
% ========================================================================
%
% BibTeX `newapa' style file for BibTeX version 0.99c, LaTeX version 2.09
% Place it in a file called newapa.bst in the BibTeX search path.  
%(Placing it in the same directory as the LaTeX document should also work.)
% Support for named citations is provided by named.sty

% This version was made by modifying the master file made by
% Oren Patashnik (PATASHNIK@SCORE.STANFORD.EDU), and the 'named' BibTeX
% style of Peter F. Patel-Schneider.

% Copyright (C) 1985, all rights reserved.
% Modifications Copyright 1989, 1990, Stephen N. Spencer
% Copying of this file is authorized only if either
% (1) you make absolutely no changes to your copy, including name, or
% (2) if you do make changes, you name it something other than
% btxbst.doc, plain.bst, unsrt.bst, alpha.bst, and abbrv.bst.
% This restriction helps ensure that all standard styles are identical.

% There are undoubtably bugs in this style.  If you make bug fixes,
% improvements, etc.  please let me know.  My e-mail address is:
%	spencer@heinlein.cgrg.ohio.state.edu
%
%   Citation format: (author-last-name, year)
%		     (author-last-name and author-last-name, year)
%		     (author-last-name {\em et al.}, year)
%                    (author-last-name)
%		     (author-last-name and author-last-name)
%		     (author-last-name {\em et al.})
%		     (year)
%
%   Reference list ordering: alphabetical by author or whatever passes
%	for author in the absence of one.
%
% This BibTeX style has support for abbreviated author lists and for
%	year-only citations.  This is done by having the citations
%	actually look like
%
%	\citeauthoryear{full-author-info}{abbrev-author-info}{year}
%
% The LaTeX style has to have the following (or similar)
%
%     \let\@internalcite\cite
%     \def\fullcite{\def\citeauthoryear##1##2##3{##1, ##3}\@internalcite}
%     \def\fullciteA{\def\citeauthoryear##1##2##3{##1}\@internalcite}
%     \def\shortcite{\def\citeauthoryear##1##2##3{##2, ##3}\@internalcite}
%     \def\shortciteA{\def\citeauthoryear##1##2##3{##2}\@internalcite}
%     \def\citeyear{\def\citeauthoryear##1##2##3{##3}\@internalcite}
%

ENTRY
  { address
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
%   month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
  }
  {}
  { label extra.label sort.label }

INTEGERS { output.state before.all mid.sentence after.sentence after.block }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.sentence :=
  #3 'after.block :=
}

STRINGS { s t u }

%%  Added this function March 9, 1990.  Used in FUNCTION format.vol.num.pages
FUNCTION {boldcize}
{ 's :=
  s "" =
    { "" }
    { "{\bf " s * "}" * }
  if$
}                           

FUNCTION {output.nonnull}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
	{ add.period$ write$
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
	    { add.period$ " " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

%% Added this entire function based on JAS98.BST  (March 10, 1990)
FUNCTION {output.nonnull1}
{ 's :=
  output.state mid.sentence =
    { ", " * write$ }
    { output.state after.block =
	{  write$                              %% This is only difference
	  newline$
	  "\newblock " write$
	}
	{ output.state before.all =
	    'write$
	    { add.period$ " " * write$ }
	  if$
	}
      if$
      mid.sentence 'output.state :=
    }
  if$
  s
}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

%% Added this function based on JAS98.BST  (March 10, 1990)
FUNCTION {output.check1}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull1               %% Only difference
  if$
}

FUNCTION {output.year.check}
{ year empty$
	 { "empty year in " cite$ * warning$ }
	 { write$
%%		" (" year * extra.label * ") " *
                ", " year * extra.label * ": " *
%%                                                      March 9, 1990
		mid.sentence 'output.state :=
	 }
  if$
}

FUNCTION {output.bibitem}
{ newline$

  "\bibitem[" write$
  label write$
  "]{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {fin.entry}
{ add.period$
  write$
  newline$
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

FUNCTION {new.sentence}
{ output.state after.block =
    'skip$
    { output.state before.all =
	'skip$
	{ after.sentence 'output.state := }
      if$
    }
  if$
}

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

FUNCTION {new.block.checka}
{ empty$
    'skip$
    'new.block
  if$
}

FUNCTION {new.block.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.block
  if$
}

FUNCTION {new.sentence.checka}
{ empty$
    'skip$
    'new.sentence
  if$
}

FUNCTION {new.sentence.checkb}
{ empty$
  swap$ empty$
  and
    'skip$
    'new.sentence
  if$
}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

FUNCTION {emphasize}
{ duplicate$ empty$
    { pop$ "" }
    { "{\em " swap$ * "}" * }
  if$
}

INTEGERS { nameptr namesleft numnames }

FUNCTION {format.names}
{ 's :=
  #1 'nameptr :=		% nameptr = 1;
  s num.names$ 'numnames :=	% numnames = num.name$(s);
  numnames 'namesleft :=
    { namesleft #0 > }

    { s nameptr "{vv~}{ll}{, jj}{, f.}" format.name$ 't :=

      nameptr #1 >
		{ namesleft #1 >
	  		{ ", " * t * }
	   		{ numnames #2 >
					{ "," * }
					'skip$
	      		if$
	      		t "others" =
						{ " et~al." * }
						{ " and " * t * }
	     	 		if$
	    	    }
	  	     if$
			 }
		    't
        if$
        nameptr #1 + 'nameptr :=			% nameptr += 1;
        namesleft #1 - 'namesleft :=			% namesleft =- 1;
    }
  while$
}

%% Added this entire function based on JAS98.BST  (March 10, 1990)
FUNCTION {format.names.reverse1}
{ 's :=
  #1 'nameptr :=		% nameptr = 1 (name pointer variable);
  s num.names$ 'numnames :=	% numnames = num.name$(s);
  numnames 'namesleft :=        % namesleft = numnames (names left variable);
    { namesleft #0 > }

%% So long as "namesleft" is positive, execute following function
    { 
%%  s nameptr "{vv~}{ll}{, jj}{, f.}" format.name$ 't :=   % Set name format

      nameptr #1 >
%%           ** Co-authors section
		{ 
     s nameptr "{f.~}{vv~}{ll}{, jj}" format.name$ 't :=   % Set name format
                        namesleft #1 >
%%                    ** Not last author
	  		{ ", " * t * }
%%                    ** Last author
	   		{ numnames #2 >
					{ "," * }
					'skip$
	      		if$
	      		t "others" =
						{ " et~al." * }
						{ " and " * t * }
	     	 		if$
	    	    }
	  	     if$
        	 }
%%               **  First author section
     { s nameptr "{vv~}{ll}{, jj}{, f.}" format.name$ 't :=   % Set name format
		    t
                    }      
        if$
        nameptr #1 + 'nameptr :=   		% nameptr += 1;
        namesleft #1 - 'namesleft :=	       	% namesleft =- 1;
    }
  while$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
%%    { author format.names }
    { author format.names.reverse1 }
%%                                              March 10, 1990
  if$
}

FUNCTION {format.key}
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names
      editor num.names$ #1 >
	{ ", editors" * }
	{ ", editor" * }
      if$
    }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }

    { title "t" change.case$ }

  if$
}

FUNCTION {n.dashify}
{ 't :=
  ""
    { t empty$ not }
    { t #1 #1 substring$ "-" =
	{ t #1 #2 substring$ "--" = not
	    { "--" *
	      t #2 global.max$ substring$ 't :=
	    }
	    {   { t #1 #1 substring$ "-" = }
		{ "-" *
		  t #2 global.max$ substring$ 't :=
		}
	      while$
	    }
	  if$
	}
	{ t #1 #1 substring$ *
	  t #2 global.max$ substring$ 't :=
	}
      if$
    }
  while$
}

FUNCTION {format.btitle}
{ title emphasize
}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {either.or.check}
{ empty$
    'pop$
    { "can't use both " swap$ * " fields in " * cite$ * warning$ }
  if$
}

FUNCTION {format.bvolume}
{ volume empty$
    { "" }
    { "volume" volume tie.or.space.connect
      series empty$
	'skip$
	{ " of " * series emphasize * }
      if$
      "volume and number" number either.or.check
    }
  if$
}

FUNCTION {format.number.series}
{ volume empty$
    { number empty$
	{ series field.or.null }
	{ output.state mid.sentence =
	    { "number" }
	    { "Number" }
	  if$
	  number tie.or.space.connect
	  series empty$
	    { "there's a number but no series in " cite$ * warning$ }
	    { " in " * series * }
	  if$
	}
      if$
    }
    { "" }
  if$
}

FUNCTION {format.edition}
{ edition empty$
    { "" }
    { output.state mid.sentence =
	{ edition "l" change.case$ " edition" * }
	{ edition "t" change.case$ " edition" * }
      if$
    }
  if$
}

INTEGERS { multiresult }

FUNCTION {multi.page.check}
{ 't :=
  #0 'multiresult :=
    { multiresult not
      t empty$ not
      and
    }
    { t #1 #1 substring$
      duplicate$ "-" =
      swap$ duplicate$ "," =
      swap$ "+" =
      or or
	{ #1 'multiresult := }
	{ t #2 global.max$ substring$ 't := }
      if$
    }
  while$
  multiresult
}

FUNCTION {format.pages}
{ pages empty$
    { "" }
    { pages multi.page.check
%%	{ "pages" pages n.dashify tie.or.space.connect }
%%	{ "page" pages tie.or.space.connect }
	{  pages n.dashify tie.or.space.connect }
	{  pages tie.or.space.connect }
      if$
    }
  if$
}

FUNCTION {format.vol.num.pages}
%%{ volume field.or.null 
{volume field.or.null boldcize                           %% March 9, 1990
  number empty$
    'skip$
    { "(" number * ")" * *
      volume empty$
	{ "there's a number but no volume in " cite$ * warning$ }
	'skip$
      if$
    }
  if$
  pages empty$
    'skip$
    { duplicate$ empty$
	{ pop$ format.pages }
%%	{ ":" * pages n.dashify * }
        { ", " * pages n.dashify * }                         %%  March 9, 1990
      if$
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages
    { type empty$
	{ "chapter" }
	{ type "l" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
	'skip$
	{ ", " * format.pages * }
      if$
    }
  if$
}

FUNCTION {format.in.ed.booktitle}
{ booktitle empty$
    { "" }
    { editor empty$
	{ "" booktitle emphasize * }
%%	{ "In " format.editors * ", " * booktitle emphasize * }
  	{ "" booktitle emphasize * ", " * format.editors * }
%%                                            March 11, 1990
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
  if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
  if$
}

FUNCTION {format.article.crossref}
{ "In"
  "\cite{" * crossref * "}" *
}

FUNCTION {format.crossref.editor}
{ editor #1 "{vv~}{ll}" format.name$
  editor num.names$ duplicate$
  #2 >
    { pop$ " et~al." * }
    { #2 <
	'skip$
	{ editor #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
	    { " et~al." * }
	    { " and " * editor #2 "{vv~}{ll}" format.name$ * }
	  if$
	}
      if$
    }
  if$
}

FUNCTION {format.book.crossref}
{ volume empty$
    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
      "In "
    }
    { "Volume" volume tie.or.space.connect
      " of " *
    }
  if$
  editor empty$
  editor field.or.null author field.or.null =
  or
    { key empty$
	{ series empty$
	    { "need editor, key, or series for " cite$ * " to crossref " *
	      crossref * warning$
	      "" *
	    }
	    { "{\em " * series * "\/}" * }
	  if$
	}
	{ key * }
      if$
    }
    { format.crossref.editor * }
  if$
  " \cite{" * crossref * "}" *
}

FUNCTION {format.incoll.inproc.crossref}
{ "In"
  " \cite{" * crossref * "}" *
}

FUNCTION {article}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% added
  output.year.check		         		% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                       %% March 10, 1990
  new.block
  crossref missing$
    { journal emphasize "journal" output.check
      format.vol.num.pages output
    }
    { format.article.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  output.year.check					% added
  new.block
%%  format.btitle "title" output.check
  format.btitle "title" output.check1                      %% March 11, 1990
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
      pages output                                   %% March 11, 1990
    }
    { new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  new.block
  note output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                          %% March 11, 1990
  new.block
  howpublished output
  address output
  pages output                                                %% March 11, 1990
  new.block
  note output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull
      crossref missing$
	{ "author and editor" editor either.or.check }
	'skip$
      if$
    }
  if$
  output.year.check						% added
  new.block
%%  format.btitle "title" output.check
  format.btitle "title" output.check1                         %% March 11, 1990
  crossref missing$
    { format.bvolume output
      new.block
      format.number.series output
      new.sentence
      publisher "publisher" output.check
      address output
%%    format.chapter.pages "chapter and pages" output.check
      pages output                                            %% March 20, 1990
    }
    { format.chapter.pages "chapter and pages" output.check
      new.block
      format.book.crossref output.nonnull
    }
  if$
  format.edition output
  new.block
  note output
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                          %% March 11, 1990
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
%%      format.chapter.pages output                           %% March 20, 1990
      new.sentence
      publisher "publisher" output.check
      address output
      format.edition output
      pages output                                            %% March 20, 1990
    }
    { format.incoll.inproc.crossref output.nonnull
      format.chapter.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                          %% March 11, 1990
  new.block
  crossref missing$
    { format.in.ed.booktitle "booktitle" output.check
      format.bvolume output
      format.number.series output
%%      new.sentence                                          %% March 17, 1990
      organization output
      publisher output
      address output
      pages output
	 }
    { format.incoll.inproc.crossref output.nonnull
      format.pages output
    }
  if$
  new.block
  note output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  format.authors output
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.btitle "title" output.check
  format.btitle "title" output.check1                      %% March 11, 1990
  organization address new.block.checkb
  organization output
  address output
  format.edition output
  pages output                                             %% March 11, 1990
  new.block
  note output
  fin.entry
}

FUNCTION {mastersthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output				% added
  output.year.check		        		% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                       %% March 11, 1990
  new.block
  "Master's thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  pages output                                             %% March 11, 1990
  new.block
  note output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  author format.key output					% added
  output.year.check						% added
  title howpublished new.block.checkb
%%  format.title output
  format.title "title" output.check1                          %% March 11, 1990
  new.block
  howpublished output
  pages output                                                %% March 11, 1990
%%  new.block                                                 %% March 17, 1990
  note output
  fin.entry
}

FUNCTION {phdthesis}
{ output.bibitem
  format.authors "author" output.check
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.btitle "title" output.check
  format.btitle "title" output.check1                         %% March 11, 1990
  new.block
  "PhD thesis" format.thesis.type output.nonnull
  school "school" output.check
  address output
  pages output                                               %% March 11, 1990
  new.block
  note output
  fin.entry
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization output }
    { format.editors output.nonnull }
  if$
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.btitle "title" output.check
  format.btitle "title" output.check1                        %% March 11, 1990
  format.bvolume output
  format.number.series output
  address output
  new.sentence
  organization output
  publisher output
  pages output                                               %% March 11, 1990
  new.block
  note output
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors "author" output.check
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                        %% March 11, 1990
  new.block
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  pages output                                              %% March 11, 1990
  new.block
  note output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output					% added
  output.year.check						% added
  new.block
%%  format.title "title" output.check
  format.title "title" output.check1                         %% March 11, 1990
  new.block
  note "note" output.check
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

INTEGERS { fullptr numfull fullsleft }

STRINGS { u1 u2 }

%FUNCTION {my.full.label}
%{ 
%  "foo" 'u :=
%  u
%}

FUNCTION {my.full.label}
{ 
% Initialize 'u1','u2','s'.
  "" 'u1 :=	   								
  "" 'u2 := 					
  's :=  					

% Initialize 'fullptr','numfull','fullsleft'.
  #1 'fullptr := 						
  s num.names$ 'numfull :=   
  numfull 'fullsleft :=		  

% enter the while loop which generates the first-citation information.
% while we have names left, 
%	 format the next name
%   if this is the next-to-last name, tack the ampersand on the end
%   else if this isn't the last name, tack the comma on the end.
%   concatenate the next name onto the first-citation string.
%   update the counters.

  { fullsleft #0 > }
  { s fullptr "{vv~}{ll}" format.name$ 'u1 :=  
	 fullsleft #2 =
	   { u1 " \& " * 'u1 := }
      { fullsleft #2 > 
	       { u1 ", " * 'u1 := }
	       'skip$
	     if$
		}
    if$
	 u2 u1 * 'u2 :=
	 fullptr #1 + 'fullptr :=  	   
	 fullsleft #1 - 'fullsleft :=  
  }
  while$

% push 'u2' onto the stack -- our first-citation information.
  u2		
}

FUNCTION {format.lab.names}
{ 's := 							
  s #1 "{vv~}{ll}" format.name$		
  s num.names$ duplicate$
  #2 >								
	 { pop$ " et~al." * }			
	 { #2 <
		  'skip$
		  { s #2 "{ff }{vv }{ll}{ jj}" format.name$ "others" =
				{ "et~al. " * }
				{ " and " * s #2 "{vv~}{ll}" format.name$ * }
			 if$
  		  }
	   if$
	 }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
	 	 { cite$ #1 #3 substring$ }
		 'key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
		  { cite$ #1 #3 substring$ }
		  'key
		if$
	 }
	 { editor format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
		  { key empty$
	    	   { cite$ #1 #3 substring$ }
	         'key
	       if$
	     }
		  { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {calc.label}
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
		  'editor.key.label
		  'author.key.label
		if$
    }
  if$
  duplicate$	

  author my.full.label  % generate the first-citation information.

  "\protect\citeauthoryear{" swap$ * "}{" * swap$ * "}{" *
  year field.or.null purify$ #-1 #4 substring$ *  
  'label :=
  year field.or.null purify$ #-1 #4 substring$ *
  sortify 'sort.label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
		  { "   " * }
	     'skip$
      if$
      s nameptr "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
		  { "et al" * }
		  { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
		  { "to sort, need author or key in " cite$ * warning$
	  	    ""
	     }
	     { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.sort}
{ editor empty$
    { key empty$
	     { "to sort, need editor or key in " cite$ * warning$
	       ""
	     }
	     { key sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { editor empty$
	     { key empty$
	         { "to sort, need author, editor, or key in " cite$ * warning$
	           ""
	         }
	         { key sortify }
	       if$
	     }
	     { editor sort.format.names }
      if$
   }
   { author sort.format.names }
  if$
}

FUNCTION {presort}
{ calc.label
  label sortify
  "    "
  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
		  'editor.sort
		  'author.sort
      if$
    }
  if$
  #1 entry.max$ substring$			% added for newapa
  'sort.label :=				% added for newapa
  sort.label					% added for newapa
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT	% sort by label, sort.label, title --- for final label calculation

STRINGS { last.label next.extra }

INTEGERS { last.extra.num }

FUNCTION {initialize.extra.label.stuff}
{ #0 int.to.chr$ 'last.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.label label =
	 { last.extra.num #1 + 'last.extra.num :=
	   last.extra.num int.to.chr$ 'extra.label :=
	 }
	 { "a" chr.to.int$ 'last.extra.num :=
	   "" 'extra.label :=
	   label 'last.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
	 'skip$
  if$
  label extra.label * "}" * 'label :=   
  extra.label 'next.extra :=
}

EXECUTE {initialize.extra.label.stuff}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {bib.sort.order}
{ sort.label
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {bib.sort.order}

SORT 		% by sort.label, year, title --- giving final bib. order.

FUNCTION {begin.bib}

{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}" write$ newline$
}

EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}
