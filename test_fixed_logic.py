
import sys
import os

# Add PKSmart to path
sys.path.append(os.path.join(os.getcwd(), "PKSmart"))

try:
    from PKSmart.app.inference import PKSmartInference
    
    print("Testing PKSmartInference...")
    engine = PKSmartInference()
    
    # Caffeine SMILES
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    result = engine.predict(smiles)
    
    print(f"SMILES: {smiles}")
    print(f"LD50 Result: {result['ld50']}")
    
    if 'value' in result['ld50'] and result['ld50']['value'] != "N/A":
        print("✅ LD50 prediction successful and transformed.")
        if result['ld50']['value'] > 0:
            print(f"✅ LD50 Value: {result['ld50']['value']} mg/kg")
        else:
            print("❌ LD50 Value is 0 or less.")
    else:
        print("❌ LD50 prediction failed (model not loaded or error).")

except Exception as e:
    print(f"❌ Test failed with error: {e}")
    import traceback
    traceback.print_exc()
