class TotalMomentumTreeData : public TObject {

 public:
  unsigned int centBinIndex;
  unsigned int particleSpeciesIndex;
  unsigned int yIndex;
  unsigned int mTm0Index;
  float avgRapidity;
  float avgmTm0;
  float avgTotalP;

  TotalMomentumTreeData() {}

  ClassDef(TotalMomentumTreeData,1);
};
