#include "Stock.h"
#include "IPricing.h"
#include "IDynamics.h"
#include "StandardMonteCarlo.h"
#include "BlackScholesDynamics.h"
#include "Barrier.h"
#include "Option.h"
#include "MertonJumpDynamics.h"



int main(){
    std::shared_ptr<MertonJumpDynamics> dynamic;
    Stock s(100, dynamic);
    DownAndOut b(ExerciseType::European , OptionType::Put, 100,100);
    


}