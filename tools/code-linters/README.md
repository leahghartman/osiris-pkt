### Collection of tools to highlight/fix/help fix various things in the code


1. _class-up-the-code_: This code uses _Fortransfigure/Wrapacious_ to go though the source
code and look for places where a polymorphic type has been improperly declared as a
Fortran ```type``` (rather then a Fortran ``class```). Example:
```

type( t_species ), intent(inout) :: this

gets changed to:

class( t_species ), intent(inout) :: this

```

By default, _class-up-the-joint_ just report on the issues it finds and DOES NOT change/write any code.
An input parameter ```--make-changes``` is needed to have the tool actually change any source code. There
are also time that are semtically unclear.. in those cases, the tool will never make changes and will,
instead, alert the problem only.
