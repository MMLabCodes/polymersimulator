class CustomDecorators:
    @staticmethod
    def call_from_method(method):
        """Decorator to set a flag indicating the method is called from another method."""
        def wrapper(instance, *args, **kwargs):
            # Check if the calling method should modify the behavior
            from_multiple = kwargs.get('from_multiple', False)  # Check for a specific keyword argument
            
            # Set the flag for indicating that the method is called from another method
            instance.called_from_another_method = True
            result = method(instance, *args, **kwargs)  # Call the actual method
            instance.called_from_another_method = False  # Reset the flag
            
            # If called from multiple, modify behavior accordingly
            if from_multiple:
                print("Behavior modified for multiple universe call.")
            return result
        return wrapper