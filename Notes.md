# Notes

- Normalization
I don't know what ype of normalization I should apply. In terms of the data, does it make sense to normalize by all samples? Are the coordinates x, y part of the same environment? What I am doing right now is to use just batch normalization per sample (i.e., not distinguishing between channels)

The problem is that, if I try any item-based normalization, I ge NaNs in the validation loss.

ACTION: try to normalize directly in the TSData object, without any ItemTransform.