#include <constraint.h>
#include <Tools/stream.h>

namespace Loci {



constraintRep::constraintRep()
{
}

constraintRep::constraintRep(const entitySet &p)
{
    constraint = p ;
}

constraintRep::~constraintRep()
{
}

void constraintRep::allocate(const entitySet &p)
{
    constraint = p ;
    dispatch_notify() ;
}

storeRep *constraintRep::new_store(const entitySet &p) const
{
    return new constraintRep(p) ;
}

store_type constraintRep::RepType() const
{
    return CONSTRAINT ;
}

const entitySet &constraintRep::domain() const {
    return constraint ;
}

ostream &constraintRep::Print(ostream &s) const {
    s << constraint << endl ;
    return s ;
}


istream &constraintRep::Input(istream &s) {
    entitySet e ;
    s >> e ;
    allocate(e) ;
    return s ;
}

constraint::constraint()
{
    setRep(new constraintType) ;
}
    
constraint::constraint(constraint &var)
{
    setRep(var.Rep()) ;
}

constraint::constraint(const entitySet &ptn)
{
    setRep(new constraintType(ptn)) ;
}

constraint::~constraint()
{
}

void constraint::notification()
{
    NPTR<constraintType> p(Rep());
    if(p!=0)
      data = p->get_constraint() ;
    warn(p==0);
}

}
