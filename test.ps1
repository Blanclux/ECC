Write-Output "*** ECC Test ***"

#cargo build --release

cargo test  --release --test ecc_test -- --nocapture
