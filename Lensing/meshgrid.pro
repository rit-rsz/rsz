function meshgrid, x,y,double=double
; duplicate the function meshgrid(x,y) in matlab
; 2011-01-30
; Ding Yuan, CFSA, Department of Physics
; Ding.Yuan@warwick.ac.uk
; 
nx=n_elements(x) & ny=n_elements(y)
if nx eq 0 then on_error,0 
if ny eq 0 then begin
y=x
ny=n_elements(y)
endif
x1=x#replicate(1,ny)
y1=y##replicate(1,nx)
; print,x
; print,y
; print,[[[x1]],[[y1]]]
result=[[[x1]],[[y1]]]
if keyword_set(double) then return,double(result) else return,result
end